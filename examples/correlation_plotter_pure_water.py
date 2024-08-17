import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from thermo import FlashPureVLS, IAPWS95Liquid, IAPWS95Gas, iapws_constants, iapws_correlations
from src.sampler.vtk_sampler import VTKSampler

folder_name = "water_ph_figures"
import os

if not os.path.exists(folder_name):
    os.makedirs(folder_name)

level = 2
file_name = 'XHP_l'+str(level)+'_iapws_modified.vtk'
brine_sampler_phz = VTKSampler(file_name)
brine_sampler_phz.conversion_factors = (1.0, 1.0, 1.0)  # (z,h,p)

# triple point of water
T_ref = 273.16
P_ref = 611.657
MW_H2O = iapws_constants.MWs[0] # [Kg/mol]
liquid = IAPWS95Liquid(T=T_ref, P=P_ref, zs=[1])
gas = IAPWS95Gas(T=T_ref, P=P_ref, zs=[1])
flasher = FlashPureVLS(iapws_constants, iapws_correlations, gas, [liquid], [])
PT = flasher.flash(P=1e7, T=573.15)

H_mass_vals = np.logspace(np.log10(1.0), np.log10(2.8), 30) # [MJ/Kg]
P_vals = np.logspace(np.log10(1.1), np.log10(20.0), 30) # [MPa]

def plot_variable(data, variable_name, title):

    z_label = {
        'Density': 'Mixture density [Kg/m3]',
        'Temperature': 'Temperature [K]',
        'Density difference': 'Mixture density [Kg/m3]',
        'Temperature difference': 'Temperature [K]'
    }

    # Example data
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    # z = vtk_T_vals

    # Create log-log plot

    # Create a 3D scatter plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c='r', marker='o')

    ax.set_xlabel('H [MJ/Kg]')
    ax.set_ylabel('P [MPa]')
    ax.set_zlabel(z_label[variable_name])
    plt.title(title)


# Density plots
PH_data = []
for H_mass in H_mass_vals:
    for P_val in P_vals:
        H_val = H_mass * MW_H2O * (1e3)  # [J/mol]
        PH = flasher.flash(P=P_val*1.0e6, H=H_val)
        PH_data.append([PH.H_mass()*1e-6,PH.P*1e-6, PH.rho_mass()])
PH_data = np.array(PH_data)

h_v = PH_data[:,0].copy()
p_v = PH_data[:,1].copy()
Rho_v = PH_data[:,2].copy()
z_NaCl = np.zeros_like(p_v)
par_points = np.array((z_NaCl, h_v, p_v)).T
brine_sampler_phz.sample_at(par_points)
vtk_rho_vals = brine_sampler_phz.sampled_could.point_data['Rho']

title = 'PH-diagram (IAPWS95)'
plot_variable(PH_data, 'Density', title)
plt.savefig(folder_name + '/density_IAPWS95.png')
plt.clf()

title = 'PH-diagram (VTKSampler)'
PH_data[:,2] = vtk_rho_vals
plot_variable(PH_data, 'Density', title)
plt.savefig(folder_name + '/density_vtk_sampler_l'+str(level)+'.png')
plt.clf()

Rho_l2_norm = np.linalg.norm(Rho_v - vtk_rho_vals)/np.linalg.norm(Rho_v)
title = 'Density difference with Rel. L2 norm = ' + str(Rho_l2_norm)
PH_data[:,2] = Rho_v - vtk_rho_vals
plot_variable(PH_data, 'Density difference', title)
plt.savefig(folder_name + '/density_difference_l'+str(level)+'.png')
plt.clf()


# Temperature plots
PH_data = []
for H_mass in H_mass_vals:
    for P_val in P_vals:
        H_val = H_mass * MW_H2O * (1e3)  # [J/mol]
        PH = flasher.flash(P=P_val*1.0e6, H=H_val)
        PH_data.append([PH.H_mass()*1e-6,PH.P*1e-6, PH.T])
PH_data = np.array(PH_data)

h_v = PH_data[:,0].copy()
p_v = PH_data[:,1].copy()
T_v = PH_data[:,2].copy()
z_NaCl = np.zeros_like(p_v)
par_points = np.array((z_NaCl, h_v, p_v)).T
brine_sampler_phz.sample_at(par_points)
vtk_T_vals = brine_sampler_phz.sampled_could.point_data['Temperature']


title = 'PH-diagram (IAPWS95)'
plot_variable(PH_data, 'Temperature', title)
plt.savefig(folder_name + '/temperature_IAPWS95.png')
plt.clf()

title = 'PH-diagram (VTKSampler)'
PH_data[:,2] = vtk_T_vals
plot_variable(PH_data, 'Temperature', title)
plt.savefig(folder_name + '/temperature_vtk_sampler_l'+str(level)+'.png')
plt.clf()

T_l2_norm = np.linalg.norm(T_v - vtk_T_vals)/np.linalg.norm(T_v)
title = 'Temperature difference with Rel. L2 norm = ' + str(T_l2_norm)
PH_data[:,2] = T_v - vtk_T_vals
plot_variable(PH_data, 'Temperature difference', title)
plt.savefig(folder_name + '/temperature_difference_l'+str(level)+'.png')
plt.clf()
