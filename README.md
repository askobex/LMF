# QuantyLF
This is a library to run ligand field calculations using Quanty, with support for multi-oxidation-state fitting.

## Features
- **Single oxidation state calculation** - Standard XAS/RIXS spectroscopy fitting
- **Multi-oxidation-state calculation** - Blend multiple oxidation states with fitted or fixed weights
- **Flexible broadening** - Support for multiple Lorentzians with negative centers and Gaussian broadening
- **Robust optimization** - Uses lmfit with Powell method for fitting

## Requirements
Python 3.10 or newer needs to be installed

## Installation
Install the package using the following command
```
$ pip install QuantyLF
```

## Example

### Single Oxidation State (Standard)

#### Prepare directory
Prepare the experimental XAS without header. The RIXS should all be combined in one file and the first lines needs to include excitation energies with preceeding 0. All values should be sperated by blanks. F.e.:

`XAS_Exp.dat`
```
629.843864705601  0.0008603675330224954
630.078187284962  0.0010035145462600522
630.260276168771  0.00034208577301302203
630.513327671302  0.0005305215326067531
...               ...
```

`RIXS_Exp.dat`
```
0 638 639.35
-5.183268710569939230e+01 0.000000e+00 0.000000e+00
-5.173604120400693063e+01 0.000000e+00 0.000000e+00
-5.163939530231446184e+01 0.000000e+00 0.0000000+00
...                       ...          ...
```


### Setup python file
An example for how the python script (`QuantyFittingLF.py`) could look like, is shown below
```py
# import package
from QuantyLF.QuantyLF import QuantyLF

# create new instance
quantyLF = QuantyLF()
# set custom quanty command
quantyLF.set_quanty_command('../Quanty_macOS', 'Darwin')

# load experimental data including excitation energies
quantyLF.load_exp_xas('XAS_Exp.dat')
quantyLF.load_exp_rixs('RIXS_Exp.dat')

# configure edge jump
# set display to True to see a plot of experimental data along edge jump (display has to be set to false for calculation)
quantyLF.config_edge_jump([[637.7, 0.14, 4], [648.2, 0.006, 8]], [600, 700,], display=False)

# see available cases
print(quantyLF.available_cases())
# load case and print parameters for fitting
quantyLF.load_case('D3h_3d', manual=True)

# Set up ion and oxidation state
quantyLF.add_par('ion', 22, from_file=False)
quantyLF.add_par('oxy', 4, from_file=False)
quantyLF.add_par('Gamma1', 0.4120250470353196, [0.4, 1])

# Crystal field contribution in D4h symmetry
quantyLF.add_par('tenDq', 0.19, [-0.2, 0.2])
quantyLF.add_par('tenDqF', 0.6815489483432551, [0.01, 1.0])

# Destruction parameter
quantyLF.add_par('Ds', 0.999, [0.1, 1])
quantyLF.add_par('Dt', 0.99, [0.1, 1])
quantyLF.add_par('DsF', 0.999, [0.1, 1])
quantyLF.add_par('DtF', 0.999, [0.1, 1])

# Multiplet contribution
# Spin orbit coupling
quantyLF.add_par('zeta_2p', 1.0196625781428472, [0.8, 1.02])
quantyLF.add_par('zeta_3d', 0.8403012992370478, [0.8, 1.02])
quantyLF.add_par('Xzeta_3d', 1.0, [0.8, 1.02])

# Slater integrals (Coulomb repulsion/exchange correlation)
quantyLF.add_par('Fdd', 0.9397729329705585, [0.8, 1.0])
quantyLF.add_par('XFdd', 0.8137253445941214, [0.8, 1.0])
quantyLF.add_par('Fpd', 0.8098173584848158, [0.8, 1.0])
quantyLF.add_par('Gpd', 0.8053014352519605, [0.8, 1.0])


# Ligand field contribution
# on-site energies (usually drops out of the equation in crystal field theory)
quantyLF.add_par('Udd', 6.543685631427877, [2.0, 7.0])
quantyLF.add_par('Upd_Udd', 4.001467895225598, [0.5, 5.0])

# Crystal field contribution of ligand site
quantyLF.add_par('tenDqL', 0.022975132006965073, [0.01, 1.0])

# Charge transfer contribution
quantyLF.add_par('Delta', 4.040660314729548, [1.0, 5.0])

# Hybridization
quantyLF.add_par('VfScale', 0.9775495515653867, [0.8, 1.0])


# # XAS and RIXS broadening
quantyLF.add_broadening('XAS', [[-3.7, 0.5], [3, 0.7], [9, 0.7]], gamma=0.5)
quantyLF.add_broadening('RIXS', [[0, 0.3], [2.7, 0.5], [8, 1]])

# Run calculation in 'RIXS' mode
quantyLF.fit('RIXS')
```

### Multi-Oxidation-State Example

Multi-oxidation-state calculation fits experimental spectra as a weighted blend of multiple oxidation states. Weights can be fitted (optimized) or fixed (user-specified).

```py
# import package
from QuantyLF.QuantyLF import QuantyLF

# create new instance
quantyLF = QuantyLF()

# load experimental data
quantyLF.load_exp_xas('XAS_Exp.dat')
quantyLF.load_exp_rixs('RIXS_Exp.dat')

# configure edge jump
quantyLF.config_edge_jump([[637.7, 0.14, 4], [648.2, 0.006, 8]], [600, 700], display=False)

# load case and set up base parameters
quantyLF.load_case('D3h_3d', manual=False)
quantyLF.add_par('ion', 22, from_file=False)

# Add multiple oxidation states
# weight=None means the weight will be fitted (optimized)
# weight=value means the weight is fixed at that value
quantyLF.add_oxidation_state(oxy=2.5, weight=None)      # Fitted weight
quantyLF.add_oxidation_state(oxy=3.0, weight=None)      # Fitted weight

# Alternative: mix fitted and fixed weights
# quantyLF.add_oxidation_state(oxy=2.0, weight=0.3)    # Fixed at 30%
# quantyLF.add_oxidation_state(oxy=2.5, weight=None)   # Fitted to ~70%

# Add other parameters (same as single oxidation state)
quantyLF.add_par('Gamma1', 0.41, [0.4, 1])
quantyLF.add_par('tenDq', 0.19, [-0.2, 0.2])
# ... add remaining parameters ...

# Broadening (supports negative centers)
quantyLF.add_broadening('XAS', [[-3.7, 0.5], [3, 0.7], [9, 0.7]], gamma=0.5)
quantyLF.add_broadening('RIXS', [[0, 0.3], [2.7, 0.5], [8, 1]])

# Run multi-oxidation-state fitting
quantyLF.fit('RIXS')

# Output files:
# - ParVals.txt: contains fitted weights (OxyWeight_2, OxyWeight_3, etc.)
# - XAS_Fit_Blended.dat: blended XAS spectrum
# - RIXS_Fit_Blended.dat: blended RIXS spectrum
```

#### Multi-Oxidation-State Workflow
1. Each oxidation state is calculated separately
2. Individual spectra are shifted and scaled to align with experimental data
3. Spectra are blended using normalized weights: `blended = Σ(weight[i] × spectrum[i])`
4. **Comparison happens on the blended spectrum only**
5. Fitted weights are optimized to minimize the difference between blended and experimental spectra

#### Weight Modes
- **Fitted weights** (`weight=None`): optimized by the fitting algorithm
- **Fixed weights** (`weight=0.3`): specified by user, not changed during fitting
- **Mixed weights**: combine fitted and fixed weights in the same calculation

**Note**: Weights are automatically normalized during fitting so they always sum to 1.0



### Running calculation using a SLURM job
Create a job file in this style:
```bash
#!/bin/bash
#SBATCH --job-name={name} # Job name
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user={mail}  # Where to send mail
#SBATCH --ntasks=1                      # Run a single task
#SBATCH --cpus-per-task=32              # Number of CPU cores per task
#SBATCH --mem=4G                   # Total memory limit
#SBATCH --time=00:30:00                 # Time limit hrs:min:sec
#SBATCH --output=pyRIXS_03G_%j.out          # Standard output and error log

# Load required modules; for example, if your program was
# compiled with Intel compiler, use the following·
# module load intel
#module load gentoo/2023
module load StdEnv/2023
module load python/3.11.5
module load scipy-stack/2023b

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index joblib
pip install -U lmfit
pip install --no-index h5py
pip install --no-index matplotlib
pip install --no-index scipy
pip install https://github.com/askobex/LMF/archive/refs/heads/main.zip


#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun python QuantyFittingLF.py #>> print_1
echo "Program finished with exit code $? at: `date`"
exit 0
```

Run calculation on computing cluster with SLURM job system
```
$ sbatch run.job
```

## Documenation

### config edge jump
```py
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
config_edge_jump(edge_jumps, x_range, y_offset=0, display=False)
```

### available_cases
```py
 """
Return the available cases for which the Quanty calculations are available

Returns
-------
List of available cases
"""
available_cases()
```

### load_case
```py
"""
Load a Quanty case

Parameters
----------
case: str
    Name of the case to load (pattern: {point_group}_{orbitals})
manual: bool, optional
    If True, the manual for the case is displayed (list of parameters available for fitting)
"""      
load_case(case, manual=False)
```

### load_custom_case
```py
"""
Load a custom Quanty case from a file

Parameters
----------
case_file: str
    Path to the Quanty case file
"""
load_custom_case(case_file)
```

### add_par
```py
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
add_par(name, init_val, interv = None, from_file = True)
```

### add_broadening
```py
"""
Add broadening to the XAS or RIXS data

Parameters
----------
type: str
    "XAS" or "RIXS"
lorenzians: list of 2-element lists
    List of lorentzian broadening parameters [center, width]
    - center: peak position on energy axis (can be negative, zero, or positive)
    - width: broadening amount (must be positive, typically 0.3-2 eV)
gamma: float, optional
    Gaussian broadening parameter (must be non-negative)

Examples
--------
# Two Lorentzians at -2 eV and +5 eV with Gaussian broadening
quantyLF.add_broadening('XAS', [[-2, 1], [5, 0.5]], gamma=0.3)

# Single Lorentzian at 0 eV
quantyLF.add_broadening('RIXS', [[0, 0.8]], gamma=0.2)
"""
add_broadening(type, lorenzians, gamma=0)
```

### add_oxidation_state
```py
"""
Add an oxidation state to the multi-oxidation-state calculation

Parameters
----------
oxy: float
    Oxidation state value (passed to Quanty calculation)
weight: float or None, optional
    Weight for this oxidation state in the blend
    - If None (default): weight is fitted (optimized)
    - If float: weight is fixed at this value (not optimized)
    
Notes
-----
- Multiple oxidation states must be added for multi-oxidation-state fitting
- Weights are automatically normalized so they sum to 1.0
- Can mix fitted (weight=None) and fixed weights in the same calculation
- Example: add_oxidation_state(oxy=2.5, weight=None) adds oxidation state 2.5 with fitted weight
- Example: add_oxidation_state(oxy=3.0, weight=0.4) adds oxidation state 3.0 with fixed weight of 0.4
"""
add_oxidation_state(oxy, weight=None)
```

### reset_oxidation_states
```py
"""
Reset all oxidation states (clear multi-oxidation-state setup)

Returns the instance to single oxidation state mode
"""
reset_oxidation_states()
```

### get_oxidation_states
```py
"""
Get list of all configured oxidation states

Returns
-------
list of tuples
    Each tuple contains (oxy_value, weight) where weight is None if fitted, or a float if fixed
"""
get_oxidation_states()

### fit
```py
"""
Fit the parameters of the model to the experimental data

Parameters
----------
mode: str
    "XAS" or "RIXS". If "XAS", only the XAS data is fitted. If "RIXS", both XAS and RIXS data is fitted

Notes
-----
Multi-Oxidation-State Fitting:
- If multiple oxidation states are added, fits blended spectrum
- Each oxidation state calculated separately
- Individual spectra shifted/scaled to align with experimental data
- Spectra blended using normalized weights
- Comparison happens on the blended spectrum only
- Output includes:
  * ParVals.txt: contains fitted weights (OxyWeight_<oxy>, etc.)
  * XAS_Fit_Blended.dat: blended XAS spectrum
  * RIXS_Fit_Blended.dat: blended RIXS spectrum (if RIXS mode)
"""
fit(type)

### load_exp_xas
```py
"""
Load the experimental XAS data from a file (no header, two columns: energy, intensity)

Parameters
----------
path: str
    Path to the file containing the experimental XAS data
"""
load_exp_xas(path)
```

### load_exp_rixs
```py
"""
Load the experimental RIXS data from a file (no header, first row: resonance energies, rest of the rows: energy, intensity)

Parameters
----------
path: str   
    Path to the file containing the experimental RIXS data
RIXS_energies: list of floats, optional
    List of resonance energies for which the RIXS data is available (if not provided, the first row of the file is used to extract the resonance energies)
"""
load_exp_rixs(path, RIXS_energies)
```

### set_quanty_command
```py
"""
Set the path to the Quanty executable (default is 'Quanty' added to path)

Parameters
----------
command: str
Path to the Quanty executable
platform: str, optional
Platform for which the path is being set (if not set the current platform is used)
"""
set_quanty_command(command, for_platform=None)
```

### set_par_value
```py
"""
Set custom path to the parameter file (default: ParVals.txt)
!Warning: Should not be changed for default cases!

Parameters
----------
par_file: str
    Path to the parameter file
"""
set_par_file(par_file)
```

```py
"""
Exports all parameters to a csv file. All scaling factors will be applied.

Parameters
----------
path: str, optional
    Path and filename to be saved to
ignore_pars: list of str, optional
    List of parameters to ignore for export
"""
export_pars(path='ExportedPars.csv', ignore_pars = ['XAS_Gamma', 'XAS_Broad', 'RIXS_Gamma', 'RIXS_Broad', 'RIXS', 'VfScale', 'Gamma1'])
```

## Develop new cases
To develop new cases, the structure of the other `.lua`-files can be used. Clone the [GitHub](https://github.com/askobex/LMF) repository and checkout the `dev` branch. Copy the file and the utils directory from [this directory](https://github.com/askobex/LMF/tree/main/src/QuantyLF/cases) into the directory of choice. Make the changes wanted. After that add the case file to the same [directory](https://github.com/askobex/LMF/tree/main/src/QuantyLF/cases). Increment the version number in the `setup.cfg` and push the `dev` branch. On GitHub create a pull request from `dev` to `main` and merge. Check if the [publishing action](https://github.com/askobex/LMF/actions) runs through.

## Multi-Oxidation-State Fitting Guide

### Overview
Multi-oxidation-state fitting allows fitting experimental XAS/RIXS spectra as a weighted blend of multiple oxidation states. This is useful when experimental data contains contributions from multiple oxidation states that coexist in the sample.

### Key Concepts

**Fitted vs Fixed Weights**
- **Fitted weight** (`weight=None`): The weight is optimized during fitting to best match experimental data
- **Fixed weight** (`weight=0.3`): The weight is set to a fixed value and not optimized

**Weight Normalization**
- Weights are automatically normalized so they always sum to exactly 1.0
- This occurs during each optimization iteration
- Allows flexible fitting without needing constraint equations

**Blending Process**
1. Calculate spectrum for each oxidation state separately
2. Shift each spectrum to align XAS peak with experimental peak
3. Scale spectrum intensity if in RIXS mode
4. Apply normalized weight to each spectrum
5. Sum weighted spectra: `blended = Σ(weight[i] × spectrum[i])`

**Comparison and Optimization**
- Comparison to experimental data happens **only on the blended spectrum**
- Individual oxidation state spectra are not compared to experimental data
- Optimizer (Powell method) minimizes the difference between blended and experimental spectra
- Fitted weights converge to physically reasonable values

### Parameter Bounds
| Parameter | Range | Notes |
|-----------|-------|-------|
| Lorentzian center | [-100, 100] eV | Position on energy axis; can be negative |
| Lorentzian width | [0.001, 20] eV | Broadening amount; must be positive |
| Gamma (Gaussian) | [0, 10] eV | Gaussian broadening; must be non-negative |
| Weights | [0.001, 1.0] | Normalized to sum to 1.0 during fitting |

### Output Files
When fitting with multiple oxidation states, the following files are generated:
- **ParVals.txt**: Contains all fitted parameters including `OxyWeight_<oxy>` for each oxidation state
- **XAS_Fit_Blended.dat**: Blended XAS spectrum (experimental vs fitted)
- **RIXS_Fit_Blended.dat**: Blended RIXS spectrum (experimental vs fitted) - only in RIXS mode
- **ParVals_oxy<N>.txt**: Individual parameter files for each oxidation state (if generated by Quanty)

### Common Issues and Solutions

**Weights not converging to expected values**
- Check that initial weights are reasonable
- Ensure oxidation state parameters are properly configured
- Verify experimental data quality
- Try adjusting optimization parameters (maxfev, tolerance)

**Broadening with negative centers**
- Negative centers are valid (represent position on energy axis)
- Example: `[[-2, 1]]` means center at -2 eV with 1 eV broadening
- Only width must be positive; center can be any value

**RIXS with multiple oxidation states**
- RIXS energies are automatically shifted by the XAS shift amount for each oxidation state
- Amplitude scaling via NNLS is applied to each state independently
- Blended RIXS spectrum includes contributions from all oxidation states

For more details, see `FINAL_IMPLEMENTATION_STATUS.md`