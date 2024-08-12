import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from thermo import FlashPureVLS, IAPWS95Liquid, IAPWS95Gas, iapws_constants, iapws_correlations
from src.sampler.vtk_sampler import VTKSampler

folder_name = "water_ph_figures"
import os

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

# triple point of water
T_ref = 273.16
P_ref = 611.657
MW_H2O = iapws_constants.MWs[0] * 1.0e-3 # [Kg/mol]
liquid = IAPWS95Liquid(T=T_ref, P=P_ref, zs=[1])
gas = IAPWS95Gas(T=T_ref, P=P_ref, zs=[1])
flasher = FlashPureVLS(iapws_constants, iapws_correlations, gas, [liquid], [])
PT = flasher.flash(P=1e7, T=573.15)

H_mass_vals = np.logspace(np.log10(100.0 * 1e3), np.log10(3000.0 * 1e3), 20) # [J/Kg]
P_vals = np.logspace(np.log10(1.0 * 1e3), np.log10(20000.0 * 1e3), 20) # [MPa]
P_vals = np.logspace(np.log10(1.1 * 1e6), np.log10(59.8 * 1e6), 20) # [MPa]
PH_data = []
T_isotherm = 300.0 + 273.15
for H_mass in H_mass_vals:
    for P_val in P_vals:
        H_val = H_mass * (18.015268*1e-3)  # [J/mol]
        PH = flasher.flash(P=P_val, H=H_val)
        PH_data.append([PH.H_mass()*1e-3,PH.P*1e-6, PH.rho_mass()])
PH_data = np.array(PH_data)


file_name = "XHP_l2_modified.vtk"
brine_sampler_phz = VTKSampler(file_name)
brine_sampler_phz.conversion_factors = (1.0, 1.0, 10.0)  # (z,h,p)

h_v = PH_data[:,0].copy()
p_v = PH_data[:,1].copy()
T_v = PH_data[:,2].copy()

z_NaCl = np.zeros_like(p_v)
par_points = np.array((z_NaCl, h_v, p_v)).T
brine_sampler_phz.sample_at(par_points)
vtk_T_vals = brine_sampler_phz.sampled_could.point_data['Rho']
diff_T_v = T_v - vtk_T_vals
diff_T_max = np.max(diff_T_v)
diff_T_min = np.min(diff_T_v)
diff_T_norm = np.linalg.norm(diff_T_v)

# Example data
x = PH_data[:,0]
y = PH_data[:,1]
z = PH_data[:,2]
# z = vtk_T_vals

# Create log-log plot

# Create a 3D scatter plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z, c='r', marker='o')

# Set axes to log scale
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_zscale('log')

ax.set_xlabel('H [KJ/Kg]')
ax.set_ylabel('P [MPa]')
ax.set_zlabel('Mixture density [Kg/m3]')
plt.title('PH-diagram (IAPWS95 z_NaCl = 0.0) ')

# # Add labels and title
# plt.xlabel('H [KJ/Kg]')
# plt.ylabel('P [MPa]')

#
# # Show grid
# plt.grid(True, which="both", ls="--")

# Display the plot
plt.show()
aka = 0
