
import numpy as np
from thermo import FlashPureVLS, IAPWS95Liquid, IAPWS95Gas, iapws_constants, iapws_correlations
from src.sampler.vtk_sampler import VTKSampler

folder_name = "water_ph_figures"
import os

from water_utils import fig_5_load_and_project_reference_data
from water_utils import fig_5_draw_and_save_comparison
from water_utils import fig_6_load_and_project_reference_data
from water_utils import fig_6_draw_and_save_comparison

if not os.path.exists(folder_name):
    os.makedirs(folder_name)


level = 2
file_name = 'XHP_l'+str(level)+'_modified_iapws.vtk'
water_sampler_phz = VTKSampler(file_name)
# z: overall composition [-]
# h: mixture/bulk enthalpy [MJ/Kg]
# p: mixture/bulk pressure [MPa]
water_sampler_phz.conversion_factors = (1.0, 1.0, 1.0)  # (z,h,p)


dx = 0.01
xc = np.arange(0.0,2.0+dx,dx)

# Moderate pressure comparison
P_proj, H_proj, T_proj, S_proj = fig_5_load_and_project_reference_data(xc)
z_proj = np.zeros_like(H_proj)
par_points = np.array((z_proj, H_proj, P_proj)).T
water_sampler_phz.sample_at(par_points)
H_vtk = water_sampler_phz.sampled_could.point_data['H']
T_vtk = water_sampler_phz.sampled_could.point_data['Temperature']
S_vtk = water_sampler_phz.sampled_could.point_data['S_l']
fig_5_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk)

# Lower pressure comparison
P_proj, H_proj, T_proj, S_proj = fig_6_load_and_project_reference_data(xc)
z_proj = np.zeros_like(H_proj)
par_points = np.array((z_proj, H_proj, P_proj)).T
water_sampler_phz.sample_at(par_points)
H_vtk = water_sampler_phz.sampled_could.point_data['H']
T_vtk = water_sampler_phz.sampled_could.point_data['Temperature']
S_vtk = water_sampler_phz.sampled_could.point_data['S_l']
fig_6_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk)
