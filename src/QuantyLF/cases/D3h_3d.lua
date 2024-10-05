-- import slater integral values
function addRelPath(dir)
    local spath_unix = debug.getinfo(1, 'S').source:sub(2):gsub("^([^/])", "./%1"):gsub("[^/]*$", "")
    print(spath_unix)
    local spath_win = debug.getinfo(1, 'S').source:sub(2):match("(.*[/\\])")
    dir = dir and (dir .. "/") or ""
    spath_unix = spath_unix .. dir
    spath_win = spath_win .. dir
    package.path = spath_unix .. "?.lua;" .. spath_unix .. "?/init.lua;" .. spath_win .. "?.lua;" .. spath_win .. "?/init.lua;" .. package.path
end

addRelPath('utils')
require "slater_integrals"

TimeStart("LF_RIXS")
-- this example calculates the resonant inelastic x-ray scattering in Transition Metal 3d we look at the
-- L-L23M45 edge, i.e. we make an excitation from 2p to 3d (L23) and decay from the
-- 3d shell back to the 2p shell (final "hole" in the 3d shell M45 edge). These spectra
-- measure d-d excitations and or magnons.

-- As magnons are not localised the final state will be dispersive. Both the energy as
-- well as the intensity will be momentum dependent. The final spectrum is thus a combination
-- of the local process as well as the non-local dispersive process.

-- Calculating a dispersive magnon in a single shot will be very hard, as one can not handle
-- large enough clusters in a full many-body calculation to capture dispersion....

-- We here take a route where we define an effective low energy operator that only allows
-- to make magnetic excitations, but has the same resonant energy and polarisation dependence
-- as the RIXS scattering. We then can use this effective operator to calculate the dispersive
-- magnons using linear spin-wave theory. (Or any other level of theory you want)
-- The local moment is alligned in the cluster due to an exchange field
HexDir = {1, 1, 1}

-- We calculate the resonant energy dependence in the following window.
Emin1 = -10
Emax1 = 20
NE1 = 120
Gamma1 = 1.0

----------------------------------------------------------------------------
----------- some functions used
----------------------------------------------------------------------------
-- calculates x^2
function sqr(x)
    return (x * x)
end

-- Read in parameters from ParVals file
function ReadParameters()

    local pars = {}
    local NSites = 0

    local citylist = {}
    for line in io.lines("ParVals.txt") do

        local iter = string.gmatch(line, "[^%s]+")
        local name = iter()
        local val = 0
        if name == "RIXS_Broad" or name == "XAS_Broad" then
            val = {tonumber(iter()), tonumber(iter())}
        else
            val = tonumber(iter())
        end
        pars[#pars + 1] = {
            name = name,
            val = val
        }
    end
    print(pars)

    return pars -- table of values, indexed by name
end

pars = ReadParameters()

NF = 26;
NB = 0;
IndexUp_2p = {0, 2, 4};
IndexDn_2p = {1, 3, 5};
IndexUp_3d = {6, 8, 10, 12, 14};
IndexDn_3d = {7, 9, 11, 13, 15};
IndexUp_Ld = {16, 18, 20, 22, 24}
IndexDn_Ld = {17, 19, 21, 23, 25}

-- Setup the ion and formal valence
ion = 0
oxy = 0

for i = 1, #pars do
    if (pars[i].name == "oxy") then
        oxy = pars[i].val
    end
    if (pars[i].name == "ion") then
        ion = pars[i].val
    end
end

-- Configure slater integrals
nd, zeta_3d, F2dd, F4dd, zeta_2p, F2pd, G1pd, G3pd, Xzeta_3d, XF2dd, XF4dd = get_slater_integrals(ion, oxy)

-- Setup the hybridization
Vb1g = 0
Va1g = 0
Vb2g = 0
Veg = 0

-- setup initial parameters for fitting values to initialize variables
Hex = 0
tenDq = 1
tenDqFs = 0
Ds = 0
Dt = 0
DsFs = 0
DtFs = 0
VfScale = 1

----overwrite init variables with ParVals values
for i = 1, #pars do
    if (pars[i].name == "Ds") then
        Ds = pars[i].val
    end
    if (pars[i].name == "Dt") then
        Dt = pars[i].val
    end
    if (pars[i].name == "DsF") then
        DsFs = pars[i].val
    end
    if (pars[i].name == "DtF") then
        DtFs = pars[i].val
    end
    if (pars[i].name == "VfScale") then
        VfScale = pars[i].val
    end
    if (pars[i].name == "Hex") then
        Hex = pars[i].val
    end
end

tenDqF = tenDqFs * tenDq
DsF = DsFs * Ds
DtF = DtFs * Dt

Va1g = 2 * Ds - 6 * Dt
Ve1g = -2 * Ds - Dt
Ve2g = -1 * Ds + 4 * Dt

Va1gF = Va1g * VfScale
Ve1gF = Ve1g * VfScale
Ve2gF = Ve2g * VfScale

OppSx_3d = NewOperator("Sx", NF, IndexUp_3d, IndexDn_3d);
OppSy_3d = NewOperator("Sy", NF, IndexUp_3d, IndexDn_3d);
OppSz_3d = NewOperator("Sz", NF, IndexUp_3d, IndexDn_3d);
OppSsqr_3d = NewOperator("Ssqr", NF, IndexUp_3d, IndexDn_3d);
OppSplus_3d = NewOperator("Splus", NF, IndexUp_3d, IndexDn_3d);
OppSmin_3d = NewOperator("Smin", NF, IndexUp_3d, IndexDn_3d);

OppLx_3d = NewOperator("Lx", NF, IndexUp_3d, IndexDn_3d);
OppLy_3d = NewOperator("Ly", NF, IndexUp_3d, IndexDn_3d);
OppLz_3d = NewOperator("Lz", NF, IndexUp_3d, IndexDn_3d);
OppLsqr_3d = NewOperator("Lsqr", NF, IndexUp_3d, IndexDn_3d);
OppLplus_3d = NewOperator("Lplus", NF, IndexUp_3d, IndexDn_3d);
OppLmin_3d = NewOperator("Lmin", NF, IndexUp_3d, IndexDn_3d);

OppJx_3d = NewOperator("Jx", NF, IndexUp_3d, IndexDn_3d);
OppJy_3d = NewOperator("Jy", NF, IndexUp_3d, IndexDn_3d);
OppJz_3d = NewOperator("Jz", NF, IndexUp_3d, IndexDn_3d);
OppJsqr_3d = NewOperator("Jsqr", NF, IndexUp_3d, IndexDn_3d);
OppJplus_3d = NewOperator("Jplus", NF, IndexUp_3d, IndexDn_3d);
OppJmin_3d = NewOperator("Jmin", NF, IndexUp_3d, IndexDn_3d);

Oppldots_3d = NewOperator("ldots", NF, IndexUp_3d, IndexDn_3d);

-- Angular momentum operators on the Ligand shell

OppSx_Ld = NewOperator("Sx", NF, IndexUp_Ld, IndexDn_Ld);
OppSy_Ld = NewOperator("Sy", NF, IndexUp_Ld, IndexDn_Ld);
OppSz_Ld = NewOperator("Sz", NF, IndexUp_Ld, IndexDn_Ld);
OppSsqr_Ld = NewOperator("Ssqr", NF, IndexUp_Ld, IndexDn_Ld);
OppSplus_Ld = NewOperator("Splus", NF, IndexUp_Ld, IndexDn_Ld);
OppSmin_Ld = NewOperator("Smin", NF, IndexUp_Ld, IndexDn_Ld);

OppLx_Ld = NewOperator("Lx", NF, IndexUp_Ld, IndexDn_Ld);
OppLy_Ld = NewOperator("Ly", NF, IndexUp_Ld, IndexDn_Ld);
OppLz_Ld = NewOperator("Lz", NF, IndexUp_Ld, IndexDn_Ld);
OppLsqr_Ld = NewOperator("Lsqr", NF, IndexUp_Ld, IndexDn_Ld);
OppLplus_Ld = NewOperator("Lplus", NF, IndexUp_Ld, IndexDn_Ld);
OppLmin_Ld = NewOperator("Lmin", NF, IndexUp_Ld, IndexDn_Ld);

OppJx_Ld = NewOperator("Jx", NF, IndexUp_Ld, IndexDn_Ld);
OppJy_Ld = NewOperator("Jy", NF, IndexUp_Ld, IndexDn_Ld);
OppJz_Ld = NewOperator("Jz", NF, IndexUp_Ld, IndexDn_Ld);
OppJsqr_Ld = NewOperator("Jsqr", NF, IndexUp_Ld, IndexDn_Ld);
OppJplus_Ld = NewOperator("Jplus", NF, IndexUp_Ld, IndexDn_Ld);
OppJmin_Ld = NewOperator("Jmin", NF, IndexUp_Ld, IndexDn_Ld);

-- total angular momentum

OppSx = OppSx_3d + OppSx_Ld
OppSy = OppSy_3d + OppSy_Ld
OppSz = OppSz_3d + OppSz_Ld
OppSsqr = OppSx * OppSx + OppSy * OppSy + OppSz * OppSz
OppLx = OppLx_3d + OppLx_Ld
OppLy = OppLy_3d + OppLy_Ld
OppLz = OppLz_3d + OppLz_Ld
OppLsqr = OppLx * OppLx + OppLy * OppLy + OppLz * OppLz
OppJx = OppJx_3d + OppJx_Ld
OppJy = OppJy_3d + OppJy_Ld
OppJz = OppJz_3d + OppJz_Ld
OppJsqr = OppJx * OppJx + OppJy * OppJy + OppJz * OppJz

-- define the coulomb operator
-- we here define the part depending on F0 seperately from the part depending on F2
-- when summing we can put in the numerical values of the slater integrals

OppF0_3d = NewOperator("U", NF, IndexUp_3d, IndexDn_3d, {1, 0, 0});
OppF2_3d = NewOperator("U", NF, IndexUp_3d, IndexDn_3d, {0, 1, 0});
OppF4_3d = NewOperator("U", NF, IndexUp_3d, IndexDn_3d, {0, 0, 1});

--- Crystal field operator for the d-shell
-- Akm = PotentialExpandedOnClm("D3h",2,{0.6, 0.0, -0.4});
-- OpptenDq_Ld = NewOperator("CF", NF, IndexUp_Ld, IndexDn_Ld, Akm)
-- Akm =PotentialExpandedOnClm("D3h",2,{ 2.0, 0.0, -7.0});
-- OppDs = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm);
-- Akm = PotentialExpandedOnClm("D3h",2,{ 4.0, 0.0, -21.0});
-- OppDt = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm);
Akm = PotentialExpandedOnClm("D3h", 2, {0.6, 0.0, -0.4});
OpptenDq_Ld = NewOperator("CF", NF, IndexUp_Ld, IndexDn_Ld, Akm)
Akm = {{2.0, 0.0, -7.0}};
OppDs = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm);
Akm = {{4.0, 0.0, -21.0}};
OppDt = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm);

-- define L-d interaction

Akm = PotentialExpandedOnClm("D3h", 2, {1, 0, 0});
OppVa1g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, Akm) +
              NewOperator("CF", NF, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, Akm)
Akm = PotentialExpandedOnClm("D3h", 2, {0, 1, 0});
OppVe1g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, Akm) +
              NewOperator("CF", NF, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, Akm)
Akm = PotentialExpandedOnClm("D3h", 2, {0, 0, 1});
OppVe2g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_Ld, IndexDn_Ld, Akm) +
              NewOperator("CF", NF, IndexUp_Ld, IndexDn_Ld, IndexUp_3d, IndexDn_3d, Akm)

OppNUp_2p = NewOperator("Number", NF, IndexUp_2p, IndexUp_2p, {1, 1, 1})
OppNDn_2p = NewOperator("Number", NF, IndexDn_2p, IndexDn_2p, {1, 1, 1})
OppN_2p = OppNUp_2p + OppNDn_2p
OppNUp_3d = NewOperator("Number", NF, IndexUp_3d, IndexUp_3d, {1, 1, 1, 1, 1})
OppNDn_3d = NewOperator("Number", NF, IndexDn_3d, IndexDn_3d, {1, 1, 1, 1, 1})
OppN_3d = OppNUp_3d + OppNDn_3d
OppNUp_Ld = NewOperator("Number", NF, IndexUp_Ld, IndexUp_Ld, {1, 1, 1, 1, 1})
OppNDn_Ld = NewOperator("Number", NF, IndexDn_Ld, IndexDn_Ld, {1, 1, 1, 1, 1})
OppN_Ld = OppNUp_Ld + OppNDn_Ld

-- Number of electrons in each of the 3d orbitals
Akm = PotentialExpandedOnClm("D3h", 2, {1, 0, 0});
OppNa1g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm)
Akm = PotentialExpandedOnClm("D3h", 2, {0, 1, 0});
OppNe1g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm)
Akm = PotentialExpandedOnClm("D3h", 2, {0, 0, 1});
OppNe2g = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, Akm)

-- In order te describe the resonance we need the interaction on the 2p shell (spin-orbit)
Oppcldots = NewOperator("ldots", NF, IndexUp_2p, IndexDn_2p);
OppUpdF0 = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {1, 0}, {0, 0});
OppUpdF2 = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 1}, {0, 0});
OppUpdG1 = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {1, 0});
OppUpdG3 = NewOperator("U", NF, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {0, 1});

-- next we define the dipole operator. The dipole operator is given as epsilon.r
-- with epsilon the polarization vector of the light and r the unit position vector
-- We can expand the position vector on (renormalized) spherical harmonics and use
-- the crystal-field operator to create the dipole operator.
-- dipole transition

t = math.sqrt(1 / 2);

Akm = {{1, -1, t}, {1, 1, -t}};
TXASx = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, Akm);
Akm = {{1, -1, t * I}, {1, 1, t * I}};
TXASy = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, Akm);
Akm = {{1, 0, 1}};
TXASz = NewOperator("CF", NF, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, Akm);

TXASr = t * (TXASx - I * TXASy);
TXASl = -t * (TXASx + I * TXASy);

TXASr.Chop()
TXASl.Chop()

TXASxdag = ConjugateTranspose(TXASx);
TXASydag = ConjugateTranspose(TXASy);
TXASzdag = ConjugateTranspose(TXASz);
TXASldag = ConjugateTranspose(TXASl);
TXASrdag = ConjugateTranspose(TXASr);

-- Setup parameters
Udd = 0
Upd = 0
tenDqL = 0
Delta = 0

-- overwrite init variables with ParVals values
for i = 1, #pars do
    if (pars[i].name == "Udd") then
        Udd = pars[i].val
    end
    if (pars[i].name == "Upd_Udd") then
        Upd = Udd + pars[i].val
    end
    if (pars[i].name == "tenDqL") then
        tenDqL = pars[i].val
    end
    if (pars[i].name == "Delta") then
        Delta = pars[i].val
    end
end

Fdd = 0.8
Fpd = 0.8
Gpd = 0.8

--- scaling with beta factor as specified in py file
for i = 1, #pars do
    if (pars[i].name == "Fdd") then
        F2dd = F2dd * pars[i].val
        F4dd = F4dd * pars[i].val
    end

    if (pars[i].name == "XFdd") then
        XF2dd = XF2dd * pars[i].val
        XF4dd = XF4dd * pars[i].val
    end

    if (pars[i].name == "Fpd") then
        F2pd = F2pd * pars[i].val
    end

    if (pars[i].name == "Gpd") then
        G1pd = G1pd * pars[i].val
        G3pd = G3pd * pars[i].val
    end

    if (pars[i].name == "zeta_3d") then
        zeta_3d = zeta_3d * pars[i].val
    end

    if (pars[i].name == "zeta_2p") then
        zeta_2p = zeta_2p * pars[i].val
    end

    if (pars[i].name == "Xzeta_3d") then
        Xzeta_3d = Xzeta_3d * pars[i].val
    end

end

F0dd = Udd + (F2dd + F4dd) * (2 / 63)
XF0dd = Udd + (XF2dd + XF4dd) * (2 / 63)
F0pd = Upd + G1pd * (1 / 15) + G3pd * (3 / 70)
Bz = 0.000001
Hz = 0.120

ed = (10 * Delta - nd * (19 + nd) * Udd / 2) / (10 + nd)
eL = nd * ((1 + nd) * Udd / 2 - Delta) / (10 + nd)

epfinal = (10 * Delta + (1 + nd) * (nd * Udd / 2 - (10 + nd) * Upd)) / (16 + nd)
edfinal = (10 * Delta - nd * (31 + nd) * Udd / 2 - 90 * Upd) / (16 + nd)
eLfinal = ((1 + nd) * (nd * Udd / 2 + 6 * Upd) - (6 + nd) * Delta) / (16 + nd)

----------------------------------------------
-- From here on it becomes different again
----------------------------------------------

-- we need two Hamiltonians, one to calculate the fundamental spectra (HEff) and one to
-- calculate the RIXS with the field in an arbitrary direction

HexDirNorm = sqrt(sqr(HexDir[1]) + sqr(HexDir[2]) + sqr(HexDir[3]))
HExchange = (Hex * HexDir[1] / HexDirNorm) * OppSx + (Hex * HexDir[2] / HexDirNorm) * OppSy +
                (Hex * HexDir[3] / HexDirNorm) * OppSz

Hamiltonian = HExchange + F0dd * OppF0_3d + F2dd * OppF2_3d + F4dd * OppF4_3d + Ds * OppDs + Dt * OppDt + zeta_3d *
                  Oppldots_3d + Bz * (2 * OppSz_3d + OppLz_3d) + Hz * OppSz_3d + tenDqL * OpptenDq_Ld + Va1g * OppVa1g +
                  Ve1g * OppVe1g + Ve2g * OppVe2g + ed * OppN_3d + eL * OppN_Ld;

XASHamiltonian = XF0dd * OppF0_3d + XF2dd * OppF2_3d + XF4dd * OppF4_3d + DsF * OppDs + DtF * OppDt + zeta_3d *
                     Oppldots_3d + Bz * (2 * OppSz_3d + OppLz_3d) + Hz * OppSz_3d + tenDqL * OpptenDq_Ld + Va1gF *
                     OppVa1g + Ve1gF * OppVe1g + Ve2gF * OppVe2g + edfinal * OppN_3d + eLfinal * OppN_Ld + epfinal *
                     OppN_2p + zeta_2p * Oppcldots + F0pd * OppUpdF0 + F2pd * OppUpdF2 + G1pd * OppUpdG1 + G3pd *
                     OppUpdG3;

-- we now can create the lowest Npsi eigenstates:
-- the calculation to the lowest 3 eigenstates.
Npsi = 3;
-- Npsi=math.fact(14)/math.fact(nd)/math.fact(14-nd)
-- in order to make sure we have a filling of 2 electrons we need to define some restrictions
res3d = "000000 1111111111 0000000000"
resL = "000000 0000000000 1111111111"
res2p = "111111 0000000000 0000000000"
StartRestrictions = {NF, NB, {res3d, nd, nd}, {res2p, 6, 6}, {resL, 10, 10}}
-- psiListEff = Eigensystem(HamiltonianEff, StartRestrictions, Npsi)
psiList = Eigensystem(Hamiltonian, StartRestrictions, Npsi, {{"restrictions", {NF, NB, {res3d, nd, nd + 1}}}})
oppList = {Hamiltonian, HExchange, OppSsqr, OppLsqr, OppJsqr, OppSz, OppLz, Oppldots_3d, OppNa1g, OppNe1g, OppNe2g,
           OppN_3d, OppN_Ld, OppN_2p};
print('\n');
print('================================================================================================\n');
print('Analysis of the initial Hamiltonian:\n');
print(' <E>  <E_ex> <S^2> <L^2>   <J^2>  <S_z> <L_z>  <l.s> N_a1g N_e1g N_e2g N_3d   N_L N_2p');
for i = 1, #psiList do
    for j = 1, #oppList do
        expectationvalue = Chop(psiList[i] * oppList[j] * psiList[i])
        io.write(string.format("%4.3f ", expectationvalue))
    end
    io.write("\n")
end
print('================================================================================================\n');
-- print( '\n');
print('=============================================\n');
print('Analysis of the Calculated orbitals energies:\n');
print("Va1g =", Va1g)
print("Ve1g =", Ve1g)
print("Ve2g =", Ve2g)
print('=============================================\n');
-- spectra XAS
doXAS = false
for i = 1, #pars do
    if (pars[i].name == "XAS") then
        doXAS = true
    end
end

if (doXAS) then

    XAS_Broad = {}
    XAS_Gamma = 0
    for i = 1, #pars do
        if (pars[i].name == "XAS_Broad") then
            XAS_Broad[#XAS_Broad + 1] = pars[i].val
        end
        if (pars[i].name == "XAS_Gamma") then
            XAS_Gamma = pars[i].val
        end
    end

    XASRestrictions = {"restrictions", {NF, NB, {res3d, nd + 1, nd + 2}}}
    -- XASSpectra = CreateSpectra(XASHamiltonian   , Tin, psiList[1],{{"Emin",-10},{"Emax",20},{"NE",2000},{"Gamma",0.1},XASRestrictions});
    Spectra_x = CreateSpectra(XASHamiltonian, TXASx, psiList[1],
        {{"Emin", -10}, {"Emax", 20}, {"NE", 2000}, {"Gamma", 0.1}, XASRestrictions});
    Spectra_z = CreateSpectra(XASHamiltonian, TXASz, psiList[1],
        {{"Emin", -10}, {"Emax", 20}, {"NE", 2000}, {"Gamma", 0.1}, XASRestrictions});
    XASSpectra = 1 / 2 * (Spectra_x + Spectra_z)
    -- workaround for case where only Gaussian broadening is used ()
    if #XAS_Broad == 0 then
        XAS_Broad = {{0, 0}}
    end
    XASSpectra.Broaden(XAS_Gamma, XAS_Broad)
    local file = assert(io.open("XAS_Calc.dat", "w"))

    SpecTables = Spectra.ToTable(XASSpectra)

    for i = 1, #SpecTables[1] do
        file:write(string.format("%14.7E ", SpecTables[1][i][1])) -- +853-peakE
        for j = 1, #SpecTables do
            file:write(string.format("%14.7E ", -SpecTables[j][i][2].Im))
        end
        file:write("\n")
    end
    file:close()
end
-- And we want to copare the effective operator RIXS calculation to the full RIXS calculation

-- spectra RIXS
doRIXS = false
RIXSEners = {}
for i = 1, #pars do
    if (pars[i].name == "RIXS") then
        doRIXS = true
        RIXSEners[#RIXSEners + 1] = pars[i].val
    end
end

RIXS_Broad = {}
for i = 1, #pars do
    if (pars[i].name == "RIXS_Broad") then
        RIXS_Broad[#RIXS_Broad + 1] = pars[i].val
    end
end

Gamma1 = 1

for i = 1, #pars do
    if (pars[i].name == "Gamma1") then
        Gamma1 = pars[i].val
    end
end

if (doRIXS) then

    RIXSRestrictions1 = {"restrictions1", {NF, NB, {res3d, nd + 1, nd + 2}}}
    RIXSRestrictions2 = {"restrictions2", {NF, NB, {resL, 9, 10}}}

    RIXSSpectra_zx = {};
    RIXSSpectra_zy = {};
    RIXSSpectra_xz = {};
    RIXSSpectra_xy = {};
    RIXSSpectra = {}
    RIXSTables = {}
    for i = 1, #RIXSEners do
        RIXSSpectra_zx[i] = CreateResonantSpectra(XASHamiltonian, Hamiltonian, TXASz, TXASxdag, psiList[1],
            {{"Emin1", RIXSEners[i]}, {"Emax1", RIXSEners[i]}, {"NE1", 1}, {"Gamma1", Gamma1}, {"Emin2", -2},
             {"Emax2", 20}, {"NE2", 2000}, {"Gamma2", 0.01}, RIXSRestrictions1, RIXSRestrictions2})
        RIXSSpectra_zy[i] = CreateResonantSpectra(XASHamiltonian, Hamiltonian, TXASz, TXASydag, psiList[1],
            {{"Emin1", RIXSEners[i]}, {"Emax1", RIXSEners[i]}, {"NE1", 1}, {"Gamma1", Gamma1}, {"Emin2", -2},
             {"Emax2", 20}, {"NE2", 2000}, {"Gamma2", 0.01}, RIXSRestrictions1, RIXSRestrictions2})
        RIXSSpectra_xz[i] = CreateResonantSpectra(XASHamiltonian, Hamiltonian, TXASx, TXASzdag, psiList[1],
            {{"Emin1", RIXSEners[i]}, {"Emax1", RIXSEners[i]}, {"NE1", 1}, {"Gamma1", Gamma1}, {"Emin2", -2},
             {"Emax2", 20}, {"NE2", 2000}, {"Gamma2", 0.01}, RIXSRestrictions1, RIXSRestrictions2})
        RIXSSpectra_xy[i] = CreateResonantSpectra(XASHamiltonian, Hamiltonian, TXASx, TXASydag, psiList[1],
            {{"Emin1", RIXSEners[i]}, {"Emax1", RIXSEners[i]}, {"NE1", 1}, {"Gamma1", Gamma1}, {"Emin2", -2},
             {"Emax2", 20}, {"NE2", 2000}, {"Gamma2", 0.01}, RIXSRestrictions1, RIXSRestrictions2})
        RIXSSpectra[i] = 0.25 * (RIXSSpectra_zx[i] + RIXSSpectra_zy[i] + RIXSSpectra_xz[i] + RIXSSpectra_xy[i])
        if #RIXS_Broad > 0 then
            RIXSSpectra[i].Broaden(0, RIXS_Broad)
        end
        RIXSTables[i] = Spectra.ToTable(RIXSSpectra[i])
    end
    -- Now write to file
    local file = assert(io.open("RIXS_Calc.dat", "w"))
    for i = 1, #RIXSTables[1][1] do -- loop over RIXS energy
        file:write(string.format("%14.7E ", RIXSTables[1][1][i][1]))
        for j = 1, #RIXSTables do -- loop over incident energy
            file:write(string.format("%14.7E ", -RIXSTables[j][1][i][2].Im))
        end
        file:write("\n")
    end
    file:close()
end

print("Finished calculating the spectra");

TimeEnd("LF_RIXS")
TimePrint()
