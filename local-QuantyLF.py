# import package
from QuantyLF.QuantyLF import QuantyLF

# create new instance
quantyLF = QuantyLF()
# set custom quanty command
quantyLF.set_quanty_command('/Users/phr542/Documents/Quanty/Quanty_macOS', 'Darwin')

# load experimental data including excitation energies
quantyLF.load_exp_xas('XAS_Exp.dat')
quantyLF.load_exp_rixs('RIXS_Exp.dat')

# configure edge jump
# set display to True to see a plot of experimental data along edge jump (display has to be set to false for calculation)
quantyLF.config_edge_jump([[637.7, 0.14, 4], [648.2, 0.006, 8]], [600, 700,], display=False)

# see available cases
# print(quantyLF.available_cases())
# load case and print parameters for fitting
quantyLF.load_case('Td_3d', manual=False)

# Set up ion and oxidation state
quantyLF.add_par('ion', 25, from_file=False)
quantyLF.add_par('oxy', 2, from_file=False)
quantyLF.add_par('Gamma1', 0.9555234791474667, [0.4, 1])

# Crystal field contribution in D4h symmetry
quantyLF.add_par('tenDq', 1.0021399464823173, [0, 2])
quantyLF.add_par('tenDqF', 0.3547998065298788, [0.01, 1.0])

# Multiplet contribution
# Spin orbit coupling
quantyLF.add_par('zeta_2p', 0.9955160828232357, [0.8, 1.02])
quantyLF.add_par('zeta_3d', 0.8454338988348356, [0.8, 1.02])
quantyLF.add_par('Xzeta_3d', 0.9241414520496289, [0.8, 1.02])

# Slater integrals (Coulomb repulsion/exchange correlation)
quantyLF.add_par('Fdd', 0.8009986743790085, [0.8, 1.0])
quantyLF.add_par('XFdd', 0.9998855686772277, [0.8, 1.0])
quantyLF.add_par('Fpd', 0.8030824105722778, [0.8, 1.0])
quantyLF.add_par('Gpd', 0.8000030853156131, [0.8, 1.0])


# Ligand field contribution
# on-site energies (usually drops out of the equation in crystal field theory)
quantyLF.add_par('Udd', 2.8041196328436917, [2.0, 7.0])
quantyLF.add_par('Upd_Udd', 1.954221327028076, [0.1, 2.0])

# Crystal field contribution of ligand site
quantyLF.add_par('tenDqL', 0.9975196109261688, [0.01, 1.0])

# Charge transfer contribution
quantyLF.add_par('Delta', 3, [1.0, 5.0], from_file=False)

# Hybridization
quantyLF.add_par('VfScale', 0.8461434962013847, [0.8, 1.0])


# # XAS and RIXS broadening
quantyLF.add_broadening('XAS', [[-3.7, 0.5], [3, 0.7], [9, 0.7]], gamma=0.5)
quantyLF.add_broadening('RIXS', [[0, 0.3], [3, 0.8], [4,0.9], [6,1]])

# Run calculation in 'RIXS' mode
# quantyLF.fit('RIXS')

quantyLF.export_pars()