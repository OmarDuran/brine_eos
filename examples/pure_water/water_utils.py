import time
import numpy as np
import pyvista as pv
import meshio
from globals import requested_fields
from thermo import FlashPureVLS, IAPWS95Liquid, IAPWS95Gas, iapws_constants, iapws_correlations
import matplotlib.pyplot as plt

def instanciate_flasher():
    # triple point of water
    T_ref = 273.16
    P_ref = 611.657
    MW_H2O = iapws_constants.MWs[0]  # [Kg/mol]
    liquid = IAPWS95Liquid(T=T_ref, P=P_ref, zs=[1])
    gas = IAPWS95Gas(T=T_ref, P=P_ref, zs=[1])
    flasher = FlashPureVLS(iapws_constants, iapws_correlations, gas, [liquid], [])
    return flasher, liquid, gas, MW_H2O

def compute_flasher_data(points):

    X, Y, Z = points.T

    H_scale = 1.0e-6
    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    # Density plots
    n_data = 13
    flasher_data = np.empty((0, n_data), float)
    for z_val, H_val, P_val in zip(X, Y, Z):
        PH = flasher.flash(P=P_val * 1e6, H=H_val * MW_H2O * 1.0e3)
        if PH.phase_count == 1:
            if PH.phase == 'L':
                data = [
                    PH.H_mass() * H_scale,
                    PH.H_mass() * H_scale,
                    PH.H_mass() * H_scale,
                    PH.rho_mass(),
                    PH.rho_mass(),
                    PH.rho_mass(),
                    1.0,
                    0.0,
                    PH.T,
                    0.0,
                    0.0,
                    PH.kinematic_viscosity(),
                    PH.kinematic_viscosity(),
                ]
                flasher_data = np.append(flasher_data, np.array([data]), axis=0)
            else:
                data = [
                    PH.H_mass() * H_scale,
                    PH.H_mass() * H_scale,
                    PH.H_mass() * H_scale,
                    PH.rho_mass(),
                    PH.rho_mass(),
                    PH.rho_mass(),
                    0.0,
                    1.0,
                    PH.T,
                    0.0,
                    0.0,
                    PH.kinematic_viscosity(),
                    PH.kinematic_viscosity(),
                ]
                flasher_data = np.append(flasher_data, np.array([data]), axis=0)
        else:
            data = [
                PH.H_mass() * H_scale,
                PH.phases[1].H_mass() * H_scale,
                PH.phases[0].H_mass() * H_scale,
                PH.rho_mass(),
                PH.phases[1].rho_mass(),
                PH.phases[0].rho_mass(),
                PH.phases[1].beta_volume,
                PH.phases[0].beta_volume,
                PH.T,
                0.0,
                0.0,
                PH.phases[1].kinematic_viscosity(),
                PH.phases[0].kinematic_viscosity(),
            ]
            flasher_data = np.append(flasher_data, np.array([data]), axis=0)
    return flasher_data


def generate_cartesian_grid(level, H_range, P_range):

    ref_data = {
        0: (100, 100),
        1: (200, 200),
        2: (300, 300),
    }

    H_vals = np.linspace(H_range[0], H_range[1], ref_data[level][0])  # [KJ/Kg]
    P_vals = np.linspace(P_range[0], P_range[1], ref_data[level][1])  # [MPa]

    x = np.array([0.0, 1.0])  # 1D array of shape (3,)
    y = H_vals  # 1D array of shape (2,)
    z = P_vals  # 1D array of shape (4,)

    # Create the meshgrid
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    # Grid dimensions
    nx, ny, nz = X.shape

    # Step 2: Create hexahedral connectivity
    # Compute the indices for the vertices of each hexahedral cell
    hexahedrons = []
    for i in range(nx - 1):
        for j in range(ny - 1):
            for k in range(nz - 1):
                # Indices of the 8 vertices of the hexahedron
                n0 = i * ny * nz + j * nz + k
                n1 = (i + 1) * ny * nz + j * nz + k
                n2 = (i + 1) * ny * nz + (j + 1) * nz + k
                n3 = i * ny * nz + (j + 1) * nz + k
                n4 = i * ny * nz + j * nz + (k + 1)
                n5 = (i + 1) * ny * nz + j * nz + (k + 1)
                n6 = (i + 1) * ny * nz + (j + 1) * nz + (k + 1)
                n7 = i * ny * nz + (j + 1) * nz + (k + 1)

                hexahedrons.append([n0, n1, n2, n3, n4, n5, n6, n7])

    # Convert to a NumPy array
    hexahedrons = np.array(hexahedrons)

    # Step 3: Create a meshio mesh
    mesh = meshio.Mesh(points, [("hexahedron", hexahedrons)])
    return mesh

def assign_point_data_with_flasher_data(mesh, flasher_data):
    points_s = mesh.points.shape
    flasher_data_s = flasher_data.shape
    assert points_s[0] == flasher_data_s[0]
    fields = {
        'H': 0,
        'H_l': 1,
        'H_v': 2,
        'Rho': 3,
        'Rho_l': 4,
        'Rho_v': 5,
        'S_l': 6,
        'S_v': 7,
        'Temperature': 8,
        'Xl': 9,
        'Xv': 10,
        'mu_l': 11,
        'mu_v': 12,
    }
    mesh.point_data['H'] = flasher_data[:, fields['H']]
    mesh.point_data['H_l'] = flasher_data[:, fields['H_l']]
    mesh.point_data['H_v'] = flasher_data[:, fields['H_v']]
    mesh.point_data['H_h'] = flasher_data[:, fields['H']]

    mesh.point_data['Rho'] = flasher_data[:, fields['Rho']]
    mesh.point_data['Rho_l'] = flasher_data[:, fields['Rho_l']]
    mesh.point_data['Rho_v'] = flasher_data[:, fields['Rho_v']]
    mesh.point_data['Rho_h'] = flasher_data[:, fields['Rho']]

    mesh.point_data['S_l'] = flasher_data[:, fields['S_l']]
    mesh.point_data['S_v'] = flasher_data[:, fields['S_v']]
    mesh.point_data['S_h'] = 0.0 * flasher_data[:, fields['S_v']]

    mesh.point_data['Temperature'] = flasher_data[:, fields['Temperature']]
    mesh.point_data['Xl'] = flasher_data[:, fields['Xl']]
    mesh.point_data['Xv'] = flasher_data[:, fields['Xv']]
    mesh.point_data['mu_l'] = flasher_data[:, fields['mu_l']]
    mesh.point_data['mu_v'] = flasher_data[:, fields['mu_v']]

def compute_gradients_and_save_vtk(level, mesh):
    # Save the mesh to a file
    file_name = 'XHP_l' + str(level) + '_modified_iapws.vtk'

    meshio.write(file_name, mesh)
    untouched_vtk_file = file_name
    print(" Computing gradients for file: ", untouched_vtk_file)
    # include gradients
    tb = time.time()
    xhp_space = pv.read(untouched_vtk_file)
    te = time.time()
    print("Loading data: Elapsed time: ", te - tb)

    tb = time.time()
    for field in requested_fields:
        if field not in xhp_space.array_names:
            raise ValueError("The vtk file does not contain the field: ", field)

    gradients = {}
    for field in requested_fields:
        # # case 1: smooths scalar field and compute gradients
        # if field in fields_to_smooth:
        #     field_data = xhp_space.point_data[field]
        #     field_sigma = 0.0
        #     if not np.isclose(np.linalg.norm(field_data),0.0):
        #         field_sigma = 0.25 * np.pi * np.std(field_data) / np.max(field_data)
        #     smooth_field_data = gaussian_filter(field_data, sigma=field_sigma)
        #     xhp_space.point_data[field] = smooth_field_data
        #
        # grad_field = xhp_space.compute_derivative(field, gradient=True)
        # gradients[field] = grad_field["gradient"]
        #
        # if field in fields_to_smooth:
        #     xhp_space.point_data[field] = field_data

        # case 2: smooths computed gradients
        # grad_field = xhp_space.compute_derivative(field, gradient=True)
        # if field in fields_to_smooth:
        #     for d in range(3):
        #         field_sigma = 0.0
        #         field_data = xhp_space.point_data[field]
        #         constant_field = np.isclose(np.linalg.norm(field_data), 0.0)
        #         constant_field_der = np.isclose(np.linalg.norm(grad_field["gradient"][:,d]), 0.0)
        #         if not (constant_field or constant_field_der):
        #             field_sigma = np.pi * np.std(field_data) / np.max(field_data)
        #         smooth_field_data = gaussian_filter(grad_field["gradient"][:,d], sigma=field_sigma)
        #         grad_field["gradient"][:,d] = grad_field["gradient"][:,d]
        #     gradients[field] = grad_field["gradient"]
        # else:
        #     gradients[field] = grad_field["gradient"]

        # case 3: no smoothing
        grad_field = xhp_space.compute_derivative(field, gradient=True)
        gradients[field] = grad_field["gradient"]

    for field in requested_fields:
        xhp_space.point_data.set_vectors(gradients[field], "grad_" + field)

    xhp_space.save(file_name, binary=True)
    te = time.time()
    print("Computing gradients on requested fields: Elapsed time: ", te - tb)

def fig_5_load_and_project_reference_data(xc):

    # doi: 10.1111/gfl.12080
    file_prefix = 'verification_pure_water/reference_data/'
    p_data = np.genfromtxt(file_prefix + 'fig_5a_pressure.csv', delimiter=',', skip_header=1)
    t_data = np.genfromtxt(file_prefix+ 'fig_5a_temperature.csv', delimiter=',', skip_header=1)
    sl_data = np.genfromtxt(file_prefix + 'fig_5b_liquid_saturation.csv', delimiter=',', skip_header=1)

    t_data[:, 1] += 273.15

    p_proj = np.interp(xc, p_data[:, 0], p_data[:, 1])
    t_proj = np.interp(xc, t_data[:, 0], t_data[:, 1])
    s_proj = np.interp(xc, sl_data[:, 0], sl_data[:, 1])

    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    def bisection(p, s, tol=1e-8, max_iter=1000):
        a = 0.0
        b = 1.0

        def func(p, s, v):
            PV = flasher.flash(P=p*1.0e6, VF=v)
            assert len(PV.betas_volume) == 2
            res = s - PV.betas_volume[0]
            return res

        if func(p, s, a) * func(p, s, b) >= 0:
            raise ValueError("f(a) and f(b) must have opposite signs")

        for _ in range(max_iter):
            c = (a + b) / 2
            if abs(func(p, s, c)) < tol or (b - a) / 2 < tol:
                return c
            elif func(p, s, c) * func(p, s, a) < 0:
                b = c
            else:
                a = c
        raise RuntimeError("Maximum iterations exceeded")

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        s_v = 1.0 - s_proj[i]
        if np.isclose(s_v, 0.0) or np.isclose(s_v, 1.0):
            PT = flasher.flash(P=pair[0]*1.0e6, T=pair[1])
            h_data.append(PT.H_mass()*1.0e-6)
        else:
            vf = bisection(pair[0], s_v)
            PT = flasher.flash(P=pair[0]*1.0e6, VF=vf)
            h_data.append(PT.H_mass()*1.0e-6)
    h_proj = np.array(h_data)
    return p_proj, h_proj, t_proj, s_proj

def fig_5_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk):
    # plot the data
    file_prefix = 'verification_pure_water/'
    figure_data = {
        'T': (file_prefix + 'figure_5_temperature_at_200_years.png', 'T - Fig. 5A P. WEIS (2014)', 'T - VTKsample + GEOMAR'),
        'S': (file_prefix + 'figure_5_liquid_saturation_at_200_years.png', 's_l - Fig. 5B P. WEIS (2014)', 's_l - VTKsample + GEOMAR'),
        'H': (file_prefix + 'figure_5_enthalpy_at_200_years.png', 'H - Fig. 5A P. WEIS (2014)', 'H - VTKsample + GEOMAR'),
    }
    fields_data = {
        'T': (T_proj,T_vtk),
        'S': (S_proj,S_vtk),
        'H': (H_proj,H_vtk),
    }

    for item in fields_data.items():
        field, data = item
        file_name, label_ref, label_vtk = figure_data[field]
        y1 = data[0]
        y2 = data[1]

        l2_norm = np.linalg.norm((data[0] - data[1])) / np.linalg.norm(data[0])

        plt.plot(xc, y1, label=label_ref)
        plt.plot(xc, y2, label=label_vtk, linestyle='--')

        plt.xlabel('Distance [Km]')
        plt.title('Relative l2_norm = ' + str(l2_norm))
        plt.legend()
        plt.savefig(file_name)
        plt.clf()

def fig_6_load_and_project_reference_data(xc):

    # doi: 10.1111/gfl.12080
    file_prefix = 'verification_pure_water/reference_data/'
    p_data = np.genfromtxt(file_prefix + 'fig_6a_pressure.csv', delimiter=',', skip_header=1)
    t_data = np.genfromtxt(file_prefix+ 'fig_6a_temperature.csv', delimiter=',', skip_header=1)
    sl_data = np.genfromtxt(file_prefix + 'fig_6b_liquid_saturation.csv', delimiter=',', skip_header=1)

    t_data[:, 1] += 273.15

    p_proj = np.interp(xc, p_data[:, 0], p_data[:, 1])
    t_proj = np.interp(xc, t_data[:, 0], t_data[:, 1])
    s_proj = np.interp(xc, sl_data[:, 0], sl_data[:, 1])

    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    def bisection(p, s, tol=1e-8, max_iter=1000):
        a = 0.0
        b = 1.0

        def func(p, s, v):
            PV = flasher.flash(P=p*1.0e6, VF=v)
            assert len(PV.betas_volume) == 2
            res = s - PV.betas_volume[0]
            return res

        if func(p, s, a) * func(p, s, b) >= 0:
            raise ValueError("f(a) and f(b) must have opposite signs")

        for _ in range(max_iter):
            c = (a + b) / 2
            if abs(func(p, s, c)) < tol or (b - a) / 2 < tol:
                return c
            elif func(p, s, c) * func(p, s, a) < 0:
                b = c
            else:
                a = c
        raise RuntimeError("Maximum iterations exceeded")

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        s_v = 1.0 - s_proj[i]
        if np.isclose(s_v, 0.0) or np.isclose(s_v, 1.0):
            PT = flasher.flash(P=pair[0]*1.0e6, T=pair[1])
            h_data.append(PT.H_mass()*1.0e-6)
        else:
            vf = bisection(pair[0], s_v)
            PT = flasher.flash(P=pair[0]*1.0e6, VF=vf)
            h_data.append(PT.H_mass()*1.0e-6)
    h_proj = np.array(h_data)
    return p_proj, h_proj, t_proj, s_proj

def fig_6_draw_and_save_comparison(xc, T_proj,T_vtk,S_proj,S_vtk,H_proj,H_vtk):
    # plot the data
    file_prefix = 'verification_pure_water/'
    figure_data = {
        'T': (file_prefix + 'figure_6_temperature_at_2000_years.png', 'T - Fig. 6A P. WEIS (2014)', 'T - VTKsample + GEOMAR'),
        'S': (file_prefix + 'figure_6_liquid_saturation_at_2000_years.png', 's_l - Fig. 6B P. WEIS (2014)', 's_l - VTKsample + GEOMAR'),
        'H': (file_prefix + 'figure_6_enthalpy_at_2000_years.png', 'H - Fig. 6A P. WEIS (2014)', 'H - VTKsample + GEOMAR'),
    }
    fields_data = {
        'T': (T_proj,T_vtk),
        'S': (S_proj,S_vtk),
        'H': (H_proj,H_vtk),
    }

    for item in fields_data.items():
        field, data = item
        file_name, label_ref, label_vtk = figure_data[field]
        y1 = data[0]
        y2 = data[1]

        l2_norm = np.linalg.norm((data[0] - data[1])) / np.linalg.norm(data[0])

        plt.plot(xc, y1, label=label_ref)
        plt.plot(xc, y2, label=label_vtk, linestyle='--')

        plt.xlabel('Distance [Km]')
        plt.title('Relative l2_norm = ' + str(l2_norm))
        plt.legend()
        plt.savefig(file_name)
        plt.clf()