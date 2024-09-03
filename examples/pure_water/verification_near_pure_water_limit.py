import time
import numpy as np
from sampler.vtk_sampler import VTKSampler

from water_utils import fig_4A_load_and_project_reference_data
from water_utils import fig_4A_draw_and_save_comparison
from water_utils import fig_4E_load_and_project_reference_data
from water_utils import fig_4E_draw_and_save_comparison
from water_utils import fig_5_load_and_project_reference_data
from water_utils import fig_5_draw_and_save_comparison
from water_utils import fig_6_load_and_project_reference_data
from water_utils import fig_6_draw_and_save_comparison

level = 2
file_name = 'XHP_l'+str(level)+'_modified_low_salt_content.vtk'
brine_sampler_phz = VTKSampler(file_name)
# z: overall composition [-]
# h: mixture/bulk enthalpy [MJ/Kg]
# p: mixture/bulk pressure [MPa]
brine_sampler_phz.conversion_factors = (1.0, 1.0e3, 10.0)  # (z,h,p)
time.sleep(7)

for k in range(10):
    dx = 0.001
    xc = np.arange(0.0,2.0+dx,dx)

    # High pressure: single phase liquid fluid
    P_proj, H_proj, T_proj = fig_4A_load_and_project_reference_data(xc)
    z_proj = np.zeros_like(H_proj)
    par_points = np.array((z_proj, H_proj, P_proj)).T
    brine_sampler_phz.sample_at(par_points)
    H_vtk = brine_sampler_phz.sampled_could.point_data['H'] * 1.0e-6
    T_vtk = brine_sampler_phz.sampled_could.point_data['Temperature']
    S_vtk = brine_sampler_phz.sampled_could.point_data['S_l']
    assert np.all(np.isclose(S_vtk, 1.0))
    fig_4A_draw_and_save_comparison(xc, T_proj,T_vtk,H_proj,H_vtk)

    # Moderate pressure: single phase steam fluid
    P_proj, H_proj, T_proj = fig_4E_load_and_project_reference_data(xc)
    z_proj = np.zeros_like(H_proj)
    par_points = np.array((z_proj, H_proj, P_proj)).T
    brine_sampler_phz.sample_at(par_points)
    H_vtk = brine_sampler_phz.sampled_could.point_data['H'] * 1.0e-6
    T_vtk = brine_sampler_phz.sampled_could.point_data['Temperature']
    S_vtk = brine_sampler_phz.sampled_could.point_data['S_l']
    assert np.all(np.isclose(S_vtk, 0.0)) # Gas states
    fig_4E_draw_and_save_comparison(xc, T_proj,T_vtk,H_proj,H_vtk)

    # Moderate pressure comparison: two-phase
    P_proj, H_proj, T_proj, S_proj = fig_5_load_and_project_reference_data(xc)
    z_proj = np.zeros_like(H_proj)
    par_points = np.array((z_proj, H_proj, P_proj)).T
    brine_sampler_phz.sample_at(par_points)
    H_vtk = brine_sampler_phz.sampled_could.point_data['H'] * 1.0e-6
    T_vtk = brine_sampler_phz.sampled_could.point_data['Temperature']
    S_vtk = brine_sampler_phz.sampled_could.point_data['S_l']
    fig_5_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk)

    # Lower pressure comparison: two-phase
    P_proj, H_proj, T_proj, S_proj = fig_6_load_and_project_reference_data(xc)
    z_proj = np.zeros_like(H_proj)
    par_points = np.array((z_proj, H_proj, P_proj)).T
    brine_sampler_phz.sample_at(par_points)
    H_vtk = brine_sampler_phz.sampled_could.point_data['H'] * 1.0e-6
    T_vtk = brine_sampler_phz.sampled_could.point_data['Temperature']
    S_vtk = brine_sampler_phz.sampled_could.point_data['S_l']
    fig_6_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk)


time.sleep(7)

brine_sampler_phz.release_memory()

time.sleep(7)
