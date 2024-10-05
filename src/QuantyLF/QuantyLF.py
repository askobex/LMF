#!/usr/bin/python
import numpy as np

from scipy.optimize import least_squares, nnls
from scipy.interpolate import interp1d
import scipy as sp
import scipy.signal as sig
import lmfit
import math
import sys
from os.path import isfile
import time
import subprocess
import platform
from .cases.utils.slater_integrals import get_slater_integrals

from multiprocessing import Process
import pickle

import matplotlib.pyplot as plt

from importlib import resources as impresources
from . import cases

class QuantyLF:

    def __init__(self):
        self.edge_jump_interp = None
        self.quanty_command = {"default": "Quanty"}
        self.platform = platform.system()
        self.par_list = []
        self.par_file = 'ParVals.txt'
        self.file_par_dict = {}
        self.__read_par_file__()

    def __read_par_file__(self):
        # check if the file exists
        if not isfile(self.par_file):
            return

        with open(self.par_file) as f:
            lines = f.readlines()
            for line in lines:
                line = line.split(maxsplit=1)
                self.file_par_dict[line[0]] = line[1].strip()

    '''
    Search through all parameters for alls matches for given name.

    Parameters
    ----------
    name: str
        Name of parameter(s) to be searched
    required: boolean, optional
        If true, exception is raised if parameter is not present at least once. Default False

    Returns
    -------
    List of parameter values that match the given name
    '''
    def __get_pars__(self, name, required = False):
        par = [par[1] for par in self.par_list if par[0] == name]
        if required and len(par) == 0:
            raise ValueError(f'No parameter found for name {name}')
        return par

    ####### Model the XAS edge jump and add it to the calculated output ###########

    def __edge_jump__(self, e, pos, jump, slope):
        return jump/np.pi * (np.arctan((e-pos)*slope) + np.pi/2)

    def __get_quanty_command__(self):
        if self.platform in self.quanty_command.keys():
            return self.quanty_command[self.platform]
        else:
            return self.quanty_command["default"]
        

    #A function to read in the edge jump in XAS as interpolated from experiment
    # def edge_jump(self, EdgeJumpObject):
    #     with open(str(EdgeJumpObject),'rb') as f:
    #         return pickle.load(f)

    def __add_edge_jump__(self, energy_scale,calcSpec):
        if self.edge_jump_interp is None:
            raise ValueError("Configure edge jump first, before running calculation with edge jump")
        edge_jump_calc = self.edge_jump_interp(energy_scale)
        calc_edge = edge_jump_calc+calcSpec
        
        return calc_edge

    ###############################################################################
    #A function which shifts a spectrum in energy by amount par[0],
    #then interpolates it to the energy spacing of the second spectrum,
    #then scales it by factor par[1], and then returns the difference
    #between the spectra

    #With edge jump fitting
    def __e_shift_res_edge_jump__(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp, add edge jump and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        calcSpecNewJump = self.__add_edge_jump__(expSpec[:,0],calcSpecNew)
        return calcSpecNewJump - expSpec[:,1]

    #Similar as above but now without edge jump fitting
    def __e_shift_res__(self, pars,calcSpec,expSpec):

        #shift calc by par, then interp and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        return calcSpecNew - expSpec[:,1]


    #Similar to above, but now returns error in derivatives - better for peak matching?
    def __e_shift_res_deriv__(self, pars,calcSpec,expSpec):

    #shift calc by par, then interp and subtract
        calcSpecNew = pars[1]*np.interp(expSpec[:,0],calcSpec[:,0]+pars[0],calcSpec[:,1])
        atanerr = np.zeros_like(calcSpecNew)
        
        for i in range(len(atanerr)):
            if(i>0):
                atanerr[i] = np.arctan2(calcSpecNew[i]-calcSpecNew[i-1],expSpec[i,0]-expSpec[i-1,0]) - np.arctan2(expSpec[i,1]-expSpec[i-1,1],expSpec[i,0]-expSpec[i-1,0])        
        return atanerr

    
    #A function which runs a Quanty calculation for a given set of parameters,
    #then compares the resulting spectrum against an experimental spectrum
    #read from file, and returns the difference between the spectra
    def __quanty_res__(self, pars,allPars,type):

        global lsiter
        
        #create a list of the names of the pars being fitted and the other pars which are not
        parNames = [x[0] for x in allPars if x[4] <= 0]
        parVals = [x[1] for x in allPars if x[4] <= 0]

        for k,v in pars.items():
            parNames.append(v.name)
            parVals.append(v.value)

        #Write the current parameters to file, so Quanty can read
        #them and do the calculation. Note for this part we make
        #sure not to write any RIXS energy parameters, because the
        #RIXS is fitted later. Here just XAS calculation is done.
        f = open("ParVals.txt","w")
        for i in range(len(parNames)):
            if(parNames[i] != "RIXS"):
                f.write(parNames[i])
                f.write(" ")
                f.write(str(parVals[i]))
                f.write("\n")   
        f.write("XAS 0\n")
        f.close()
        

        #Print out current values of parameters, so they can be tracked
        print("=====================")
        for i in range(len(parNames)):
            print(parNames[i],parVals[i])
        print("=====================")


        #run Quanty - it will read the parameters from file, calculate 
        #the XAS, and write the new spectrum/spectra to file XAS_Calc.dat
        # subprocess.call([self.__get_quanty_command__(), self.quanty_file],stdout=subprocess.DEVNULL)
        subprocess.call([self.__get_quanty_command__(), self.quanty_file])#,stdout=None)

        
        #load spectra and experiment to compare
        calcSpec = np.loadtxt("XAS_Calc.dat")        
        
        #find the energy of largest peak in calculated and
        #experiment, to give a first rough shift for aligning energy
        calcPeak = calcSpec[calcSpec[:,1].argmax(),0]
        expPeak = self.expXAS[self.expXAS[:,1].argmax(),0]
        
        #parameters for the least squares function to fit the shift
        #and amplitude of the calculated XAS to compare to experiment
        dE = np.array([expPeak-calcPeak, 1]) #need to give it a close guess for energy shift
        amp = np.array([1])
        lowlim = np.array([-1e5,0])
        highlim = np.array([1e5,np.inf])
        
        #perform a non-linear least squares fit of the energy shift and amplitude of the 
        #calculated XAS, in order to compare it to experiment
        #Two versions are below - one using regular difference for least squares, and one using derivatives
        
        res_fn = self.__e_shift_res_edge_jump__ if self.edge_jump_interp is not None else self.__e_shift_res__

        res = least_squares(res_fn,dE,bounds=(lowlim,highlim),max_nfev=200,args=(calcSpec,self.expXAS),verbose=0)#

        #Get the difference between calculated and experiment
        diff = res.fun
        diff_XAS = np.true_divide(diff,len(self.expXAS))
        
        #Write the XAS to file (this is on the calculated energy grid)
        calcSpec[:,0] = calcSpec[:,0] + res.x[0]
        calcSpec[:,1] = calcSpec[:,1] * res.x[1]
        
        if self.edge_jump_interp is not None:
            calcSpec[:,1] = self.__add_edge_jump__(calcSpec[:,0],calcSpec[:,1])
        
        np.savetxt("XAS_Fit.dat",calcSpec)


        #If requested by the user, now do a fit of parameters based
        #on RIXS spectra
        if(type == "RIXS"):

            #rewrite the par file now with RIXS energies included, shifted appropriately
            f = open("ParVals.txt","w")
            RIXSEner = []
            for i in range(len(parNames)):
                f.write(parNames[i])
                f.write(" ")
                if(parNames[i]=="RIXS"):
                    print(parNames[i],parVals[i])
                    f.write(str(parVals[i]-res.x[0])) #shifted resonance energy according to XAS shift determined above
                    RIXSEner.append(parVals[i])
                else:
                    f.write(str(parVals[i])) #other parameter (non RIXS)      
                f.write("\n")
            f.close()    
            
            #call Quanty to do the RIXS calculation with the current set of parameters
            subprocess.call([self.__get_quanty_command__(), self.quanty_file])#,stdout=subprocess.DEVNULL)

            #load the calculated RIXS spectra
            calcRIXS = np.loadtxt("RIXS_Calc.dat")
            
            
            #fit scaling of the RIXS spectra to get best agreement
            #with experiment (linear least squares fit, single scaling for all)
            #to do this, concatenate all the RIXS spectra
            calcRIXS2 = np.copy(self.expRIXS)
            for i in range(len(calcRIXS2[0,:])-1):
                calcRIXS2[:,i+1] = np.interp(self.expRIXS[:,0],calcRIXS[:,0],calcRIXS[:,i+1])
            calcRIXSCat = np.copy(calcRIXS2[:,1])
            expRIXSCat = np.copy(self.expRIXS[:,1])
            for i in range(len(calcRIXS[0,:])-2):
                calcRIXSCat = np.hstack((calcRIXSCat,calcRIXS2[:,i+2]))
                expRIXSCat = np.hstack((expRIXSCat,self.expRIXS[:,i+2]))
            
            #do the linear least squares fit
            amp,res = nnls(np.array(calcRIXSCat[:,None]),expRIXSCat)

            print(amp, res)
            
            #Apply the scaling factor to the calculated RIXS
            #(Both the interpolated and the original calculated)
            for i in range(len(calcRIXS2[0,:])-1):
                calcRIXS2[:,i+1] = amp[0] * calcRIXS2[:,i+1]
                calcRIXS[:,i+1] = amp[0] * calcRIXS[:,i+1]

            #save RIXS, return concatenated differences
            np.savetxt("RIXS_Fit.dat",calcRIXS)
            np.savetxt("RIXS_Exp_Trimmed.dat",self.expRIXS)
            diff = calcRIXS2[:,1] - self.expRIXS[:,1]
            for i in range(len(calcRIXS2[0,:])-2):
                diff = np.hstack((diff,calcRIXS2[:,i+2] - self.expRIXS[:,i+2]))
            print("Chi2: ",np.dot(diff,diff))
            sys.stdout.flush()
            
            diff_RIXS = np.true_divide(abs(calcRIXS2[:,1] - self.expRIXS[:,1]),len(self.expRIXS))
            counter = 1
            for i in range(len(calcRIXS2[0,:])-2):
                diff_RIXS = np.hstack((diff_RIXS,np.true_divide(abs(calcRIXS2[:,i+2] - self.expRIXS[:,i+2]),len(self.expRIXS))))
                counter += 1
            
            diff_RIXS2 = np.dot(diff_RIXS,diff_RIXS)/counter
            diff_XAS2 = np.dot(diff_XAS,diff_XAS)
            
            if lsiter == 1:
                global initialXAS
                global initialRIXS
                
                initialXAS = diff_XAS2
                initialRIXS = diff_RIXS2
                                
            diff_weighted = 0.5*(diff_XAS2/initialXAS) + 0.5*(diff_RIXS2/initialRIXS)
                                    
            print("Diff XAS",diff_XAS2)
            print("Diff RIXS",diff_RIXS2)
            print("Diff weighted",diff_weighted)
            lsiter +=1
                                    
            return diff_weighted

            
        #XAS fitting, so return XAS error
        print("Chi2: ",np.dot(diff,diff))
        sys.stdout.flush()
        
        return  diff #should be able to return something like res.error? 

    
    def __fit_pars__(self,type):
        params = lmfit.Parameters()
        for i in self.par_list:
            if i[4] == 1:
                params.add(i[0],value=i[1],min=i[2],max=i[3])

    #QuantyRes(pars,type)
        global lsiter
        lsiter = 1

        minimizer = lmfit.Minimizer(self.__quanty_res__,params,fcn_args=(self.par_list,type))
        res = minimizer.minimize(method='powell',params=params,options={'xtol': 1e-12})

        print("Final values: ")
        print(res.params)



    """
    Configure the edge jump for the XAS calculation. The edge jump is modelled as a arctan function.

    Parameters
    ----------
    edge_jumps: list of 3-element lists
        List of edge jumps. Each element in the list is a list of 3 elements: [position, jump, slope]
    x_range: list of floats
        The range of x values to model the edge jump over
    y_offset: float, optional
        The y offset of the edge jump
    display: bool, optional
        If True, the edge jump will be displayed along with the experimental XASs
    """
    def config_edge_jump(self, edge_jumps, x_range, y_offset=0, display=False):
        x_range = np.linspace(x_range[0], x_range[1], math.ceil((max(x_range)-min(x_range))/0.1))
        edge_jump_y = [y_offset] * len(x_range)
        for i in range(len(edge_jumps)):
            cur_jump = edge_jumps[i]
            edge_jump_y += self.__edge_jump__(x_range, cur_jump[0], cur_jump[1], cur_jump[2])
        
        self.edge_jump_interp = interp1d(x_range, edge_jump_y, kind='cubic', fill_value="extrapolate")

        if display:
            print("Displaying edge jump, this option has to be disabled for batch runs!")
            if self.expXAS is None:
                print("To display edge jump with experimental XAS, load the experimental XAS first")
            else:
                plt.plot(self.expXAS[:,0], self.expXAS[:,1])
            plt.plot(x_range, edge_jump_y)
            plt.show()

    """
    Return the available cases for which the Quanty calculations are available

    Returns
    -------
    List of available cases
    """
    def available_cases(self):
        base_path = impresources.files(cases)
        available_cases = []
        for file in base_path.iterdir():
            if file.suffix == '.lua':
                available_cases.append(file.stem)

        return available_cases
    

    """
    Load a Quanty case

    Parameters
    ----------
    case: str
        Name of the case to load (pattern: {point_group}_{orbitals})
    manual: bool, optional
        If True, the manual for the case is displayed (list of parameters available for fitting)
    """    
    def load_case(self, case, manual=False):
        base_path = impresources.files(cases)
        case_path = base_path / f"{case}.lua"
        if not case_path.exists():
            raise ValueError(f"Case {case} not found")
        self.quanty_file = case_path

        if manual:
            manual_path = base_path / f"{case}.txt"
            if manual_path.exists():
                print(f"Case {case} loaded")
                print(f'The following parameters are available for fitting of {case}')
                with open(manual_path) as f:
                    print(f.read())
            else:
                print(f"Case {case} loaded, but no manual available")
            


    """
    Load a custom Quanty case from a file

    Parameters
    ----------
    case_file: str
        Path to the Quanty case file
    """
    def load_custom_case(self, case_file):
        if not isfile(case_file):
            raise ValueError(f"File {case_file} not found")
        self.quanty_file = case_file



    """
    Add a parameter to the model

    Parameters
    ----------
    name: str
        Name of the parameter
    init_val: float
        Initial value of the parameter
    interv: list of 2 floats, optional
        List of two values, lower and upper bounds of the parameter
    from_file: bool, optional
        If True (default True), the parameter is read from a file (file value overrides init_val)
    """
    def add_par(self, name, init_val, interv = None, from_file = True):
        if interv and (interv[1] - interv[0]) < 0:
            raise ValueError("Upper bound of parameter should be greater than lower bound")
        if interv and (init_val < interv[0] or init_val > interv[1]):
            raise ValueError("Initial value of parameter should be within the bounds")

        low = 0
        high = 0
        if interv:
            low = interv[0]
            high = interv[1]
            if name in self.file_par_dict.keys():
                init_val = float(self.file_par_dict[name]) if from_file else init_val
        else:
            if name in self.file_par_dict.keys():
                init_val = self.file_par_dict[name] if from_file else init_val
        self.par_list.append([name, init_val, low, high, 1 if interv else 0])


    """
    Add broadening to the XAS or RIXS data

    Parameters
    ----------
    type: str
        "XAS" or "RIXS"
    lorenzians: list of 2-element lists
        List of lorenzian broadening parameters (center, width)
    gamma: float, optional
        Guassian broadening parameter
    """
    def add_broadening(self, type, lorenzians, gamma=0):
        self.add_par(type + "_Gamma", gamma)
        for lorenzian in lorenzians:
            val = f'{lorenzian[0]} {lorenzian[1]}'
            self.add_par(type + "_Broad", val, from_file=False)

    """
    Fit the parameters of the model to the experimental data

    Parameters
    ----------
    mode: str
        "XAS" or "RIXS". If "XAS", only the XAS data is fitted. If "RIXS", both XAS and RIXS data is fitted
    """
    def fit(self, type):        
        self.__fit_pars__(type)

    """
    Load the experimental XAS data from a file (no header, two columns: energy, intensity)

    Parameters
    ----------
    path: str
        Path to the file containing the experimental XAS data
    """
    def load_exp_xas(self, path):
        self.expXAS = np.loadtxt(path)


    """
    Load the experimental RIXS data from a file (no header, first row: resonance energies, rest of the rows: energy, intensity)

    Parameters
    ----------
    path: str   
        Path to the file containing the experimental RIXS data
    RIXS_energies: list of floats, optional
        List of resonance energies for which the RIXS data is available (if not provided, the first row of the file is used to extract the resonance energies)
    """
    def load_exp_rixs(self, path, RIXS_energies=None):        
        #load exp RIXS. For experimental, the first row are the resonance energies
        expRIXS = np.loadtxt(path)
        energies = expRIXS[0,:]
        if energies[0] != 0:
            raise ValueError("First value of the first row should be 0, as it is the energy axis for the RIXS data")        
        if RIXS_energies is None:
            RIXS_energies = energies[1:]
        else:
            # check if number rixs energies is equal to number of columns in expRIXS
            if len(RIXS_energies) != len(expRIXS[1,:])-1:
                raise ValueError("Number of RIXS energies does not match the number of columns in the experimental RIXS data")
            # check if the rixs energies match energies in the first row of the expRIXS
            for i in range(len(RIXS_energies)):
                if RIXS_energies[i] != energies[i+1]:
                    raise ValueError(f"Resonant energy {RIXS_energies[i]} not found in the experimental RIXS data. First row should contain the resonance energies")

        for RIXS_energy in RIXS_energies:
            self.add_par("RIXS",RIXS_energy,from_file=False)
        
        #trim the exp RIXS just to have the columns with resonant energies we are calculating
        indices = [0] #0th column is energy, which we will keep
        for i in range(len(RIXS_energies)):
            for j in range(len(expRIXS[0,:])):
                if(abs(RIXS_energies[i]-expRIXS[0,j]) < 0.1):
                    indices.append(j)
            
        self.expRIXS = expRIXS[1:,indices] #this removes the first row as well, which had the energy values (which we no longer need)


    """
    Set the path to the Quanty executable (default is 'Quanty' added to path)

    Parameters
    ----------
    command: str
        Path to the Quanty executable
    platform: str, optional
        Platform for which the path is being set (if not set the current platform is used)
    """
    def set_quanty_command(self, command, for_platform=None):
        if for_platform is None:
            for_platform = platform.system()
        self.quanty_command[for_platform] = command


    """
    Set custom path to the parameter file (default: ParVals.txt)
    !Warning: Should not be changed for default cases!

    Parameters
    ----------
    par_file: str
        Path to the parameter file
    """
    def set_par_file(self, par_file):
        self.par_file = par_file



    """
    Exports all parameters to a csv file. All scaling factors will be applied.

    Parameters
    ----------
    path: str, optional
        Path and filename to be saved to, default './ExportedPars.csv'
    ignore_pars: list of str, optional
        List of parameters to ignore for export, default set to ['XAS_Gamma', 'XAS_Broad', 'RIXS_Gamma', 'RIXS_Broad', 'RIXS', 'VfScale', 'Gamma1']
    """
    def export_pars(self, path='ExportedPars.csv', ignore_pars = ['XAS_Gamma', 'XAS_Broad', 'RIXS_Gamma', 'RIXS_Broad', 'RIXS', 'VfScale', 'Gamma1']):
        ion = self.__get_pars__('ion', required=True)[0]
        oxy = self.__get_pars__('oxy', required=True)[0]

        zeta_3d, F2dd, F4dd, zeta_2p, F2pd, G1pd, G3pd, Xzeta_3d, XF2dd, XF4dd = get_slater_integrals(ion, oxy)

        export_par_list = []
        for par in self.par_list:
            name = par[0]
            val = par[1]
            if name in ignore_pars:
                continue

            # apply scaling factors
            match name:
                case 'tenDqF':
                    tenDq = self.__get_pars__('tenDq', required=True)[0]
                    export_par_list.append([name, float(val) * float(tenDq)])
                case 'zeta_2p':
                    export_par_list.append([name, float(val) * float(zeta_2p)])
                case 'zeta_3d':
                    export_par_list.append([name, float(val) * float(zeta_3d)])
                case 'Xzeta_3d':
                    export_par_list.append([name, float(val) * float(Xzeta_3d)])
                case 'Fdd':
                    export_par_list.append(['F2dd', float(val) * float(F2dd)])
                    export_par_list.append(['F4dd', float(val) * float(F4dd)])
                case 'XFdd':
                    export_par_list.append(['XF2dd', float(val) * float(XF2dd)])
                    export_par_list.append(['XF4dd', float(val) * float(XF4dd)])               
                case 'Fpd':
                    export_par_list.append(['F2pd', float(val) * float(F2pd)])
                case 'Gpd':
                    export_par_list.append(['G1pd', float(val) * float(G1pd)])
                    export_par_list.append(['G3pd', float(val) * float(G3pd)])
                case _:
                    export_par_list.append([name, val])

        export_par_list = np.asarray(export_par_list)
        np.savetxt('ExportedPars.csv', export_par_list, delimiter=',', fmt="%s")
        