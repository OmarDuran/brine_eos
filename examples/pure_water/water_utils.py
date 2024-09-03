import time
import numpy as np
import pyvista as pv
import gmsh
import sys
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
    counter = 0
    for z_val, H_val, P_val in zip(X, Y, Z):
        print('ip: ', counter)
        print('(H_val, P_val): ', (H_val, P_val))
        liquid_q = False
        vapor_q = False
        try:
            PH = flasher.flash(P=P_val * 1e6, H=H_val * MW_H2O * 1.0e3)
        except:
            H_spec = H_val * MW_H2O * 1.0e3
            HV_liquid = flasher.flash(P=P_val * 1e6, VF=0.0)
            HV_vapor = flasher.flash(P=P_val * 1e6, VF=1.0)
            if np.isclose(H_spec,HV_liquid.H()):
                PH = HV_liquid
                liquid_q = True
            elif np.isclose(H_spec,HV_vapor.H()):
                PH = HV_vapor
                vapor_q = True
            else:
                raise ValueError('compute_flasher_data:: Unable to compute thermodynamic state.')
        counter +=1
        two_phase_data = len(PH.betas_mass) == 2
        if two_phase_data:
            if np.isclose(PH.betas_mass[0], 0.0):
                liquid_q = True
            else:
                vapor_q = True
        if len(PH.betas_mass) == 1:
            if PH.phase == 'L' or liquid_q:
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
                    PH.mu(),
                    PH.mu(),
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
                    PH.mu(),
                    PH.mu(),
                ]
                flasher_data = np.append(flasher_data, np.array([data]), axis=0)
        elif two_phase_data and (liquid_q or vapor_q):
            if liquid_q:
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
                    PH.mu(),
                    PH.mu(),
                ]
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
                    PH.mu(),
                    PH.mu(),
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
                PH.phases[1].mu(),
                PH.phases[0].mu(),
            ]
            flasher_data = np.append(flasher_data, np.array([data]), axis=0)
    return flasher_data


def generate_cartesian_grid(level, H_range, P_range):

    ref_data = {
        0: (100, 100),
        1: (200, 200),
        2: (400, 400),
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

def generate_simplex_grid_pure_water(level, H_range, P_range):

    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    lc_scale = 1.0
    ref_data = {
        0: (5, 200),
        1: (100, 400),
        2: (200, 800),
    }

    H_vals = np.linspace(H_range[0], H_range[1], ref_data[level][0])  # [KJ/Kg]
    P_vals = np.linspace(P_range[0], P_range[1], ref_data[level][1])  # [MPa]

    H_scale = 1.0e-6
    H_vapor = []
    H_liquid = []
    for P_val in P_vals:
        PL = flasher.flash(P=P_val * 1.0e6, VF=0.0, zs=[1.0])
        PV = flasher.flash(P=P_val * 1.0e6, VF=1.0, zs=[1.0])
        H_liquid.append(PL.H_mass() * H_scale)
        H_vapor.append(PV.H_mass() * H_scale)

    # Define points for the polylines inside the rectangular domain
    curve_liquid = np.vstack([ np.zeros_like(P_vals), H_liquid, P_vals]).T
    curve_vapor = np.vstack([ np.zeros_like(P_vals), H_vapor, P_vals]).T

    # Initialize Gmsh
    gmsh.initialize()
    gmsh.model.add("IAPWS simplex polyline Enthalpy-Pressure space")

    lc = np.abs((H_range[1] - H_range[0]) / ref_data[level][0])

    # Define the points for the rectangular domain
    p1 = gmsh.model.occ.addPoint(0, H_range[0], P_range[0], lc)
    p2 = gmsh.model.occ.addPoint(0, H_range[1], P_range[0], lc)
    p3 = gmsh.model.occ.addPoint(0, H_range[1], P_range[1], lc)
    p4 = gmsh.model.occ.addPoint(0, H_range[0], P_range[1], lc)

    # points from the phase boundary
    p5 = gmsh.model.occ.addPoint(curve_liquid[0][0], curve_liquid[0][1],
                                 curve_liquid[0][2], lc)
    p6 = gmsh.model.occ.addPoint(curve_vapor[0][0], curve_vapor[0][1],
                                 curve_vapor[0][2], lc)
    p7 = gmsh.model.occ.addPoint(curve_liquid[-1][0], curve_liquid[-1][1],
                                 curve_liquid[-1][2], lc)
    p8 = gmsh.model.occ.addPoint(curve_vapor[-1][0], curve_vapor[-1][1],
                                 curve_vapor[-1][2], lc)


    # Define the lines for the rectangular domain including segments of intersections
    l1 = gmsh.model.occ.addLine(p1, p5)
    l2 = gmsh.model.occ.addLine(p5, p6)
    l3 = gmsh.model.occ.addLine(p6, p2)
    l4 = gmsh.model.occ.addLine(p2, p3)
    l5 = gmsh.model.occ.addLine(p3, p8)
    l6 = gmsh.model.occ.addLine(p8, p7)
    l7 = gmsh.model.occ.addLine(p7, p4)
    l8 = gmsh.model.occ.addLine(p4, p1)

    # Create a curve loop and a surface for the rectangular domain
    rect_loop = gmsh.model.occ.addCurveLoop([l1, l2, l3, l4, l5, l6, l7, l8])
    rect_surface = gmsh.model.occ.addPlaneSurface([rect_loop])

    # Add points and curves for liquid
    curve_liquid_ids = []
    for x, y, z in curve_liquid:
        point_id = gmsh.model.occ.addPoint(x, y, z, lc_scale * lc)
        curve_liquid_ids.append(point_id)

    curve_liquid_ids[0] = p5
    curve_liquid_ids[-1] = p7
    # Create lines between consecutive points to form the polyline
    curve_liquid_lines = []
    for i in range(len(curve_liquid_ids) - 1):
        line_id = gmsh.model.occ.addLine(curve_liquid_ids[i], curve_liquid_ids[i + 1])
        curve_liquid_lines.append(line_id)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.embed(1, curve_liquid_lines, 2, rect_surface)

    # Add points and curves for liquid
    curve_vapor_ids = []
    for x, y, z in curve_vapor:
        point_id = gmsh.model.occ.addPoint(x, y, z, lc_scale * lc)
        curve_vapor_ids.append(point_id)
    curve_vapor_ids[0] = p6
    curve_vapor_ids[-1] = p8
    # Create lines between consecutive points to form the polyline
    curve_vapor_lines = []
    for i in range(len(curve_vapor_ids) - 1):
        line_id = gmsh.model.occ.addLine(curve_vapor_ids[i], curve_vapor_ids[i + 1])
        curve_vapor_lines.append(line_id)
    gmsh.model.occ.synchronize()
    gmsh.model.mesh.embed(1, curve_vapor_lines, 2, rect_surface)
    gmsh.model.occ.synchronize()

    # gmsh.model.addPhysicalGroup(2, [rect_surface], name="Rectangular Surface")
    # gmsh.model.addPhysicalGroup(1, curve_liquid_ids, name="Liquid boundary")
    # gmsh.model.addPhysicalGroup(1, curve_vapor_ids, name="Vapor boundary")
    gmsh.model.mesh.generate(2)
    gmsh.write("polylines_in_rectangle.msh")
    gmsh.finalize()

    mesh = meshio.read("polylines_in_rectangle.msh")
    mesh.cell_sets.pop('gmsh:bounding_entities', None)
    mesh.cell_data.pop('gmsh:geometrical', None)
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

    mesh.point_data.pop('gmsh:dim_tags', None)

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
        # case 3: no smoothing
        grad_field = xhp_space.compute_derivative(field, gradient=True)
        gradients[field] = grad_field["gradient"]

    for field in requested_fields:
        xhp_space.point_data.set_vectors(gradients[field], "grad_" + field)

    xhp_space.save(file_name, binary=True)
    te = time.time()
    print("Computing gradients on requested fields: Elapsed time: ", te - tb)

def fig_4A_load_and_project_reference_data(xc):

    # doi: 10.1111/gfl.12080
    file_prefix = 'verification_pure_water/reference_data/'
    p_data = np.genfromtxt(file_prefix + 'fig_4a_pressure.csv', delimiter=',', skip_header=1)
    t_data = np.genfromtxt(file_prefix+ 'fig_4a_temperature.csv', delimiter=',', skip_header=1)

    t_data[:, 1] += 273.15
    p_proj = np.interp(xc, p_data[:, 0], p_data[:, 1])
    t_proj = np.interp(xc, t_data[:, 0], t_data[:, 1])

    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        PT = flasher.flash(P=pair[0] * 1.0e6, T=pair[1])
        h_data.append(PT.H_mass() * 1.0e-6)
    h_proj = np.array(h_data)
    return p_proj, h_proj, t_proj

def fig_4A_draw_and_save_comparison(xc, T_proj,T_vtk,H_proj,H_vtk):
    # plot the data
    file_prefix = 'verification_pure_water/'
    figure_data = {
        'T': (file_prefix + 'figure_4a_temperature_at_250_years.png', 'T - Fig. AA P. WEIS (2014)', 'T - VTKsample + GEOMAR'),
        'H': (file_prefix + 'figure_4a_enthalpy_at_250_years.png', 'H - Fig. AA P. WEIS (2014)', 'H - VTKsample + GEOMAR'),
    }
    fields_data = {
        'T': (T_proj,T_vtk),
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

def fig_4E_load_and_project_reference_data(xc):

    # doi: 10.1111/gfl.12080
    file_prefix = 'verification_pure_water/reference_data/'
    p_data = np.genfromtxt(file_prefix + 'fig_4e_pressure.csv', delimiter=',', skip_header=1)
    t_data = np.genfromtxt(file_prefix+ 'fig_4e_temperature.csv', delimiter=',', skip_header=1)

    t_data[:, 1] += 273.15
    p_proj = np.interp(xc, p_data[:, 0], p_data[:, 1])
    t_proj = np.interp(xc, t_data[:, 0], t_data[:, 1])

    flasher, liquid, gas, MW_H2O = instanciate_flasher()

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        PT = flasher.flash(P=pair[0] * 1.0e6, T=pair[1])
        h_data.append(PT.H_mass() * 1.0e-6)
    h_proj = np.array(h_data)
    return p_proj, h_proj, t_proj

def fig_4E_draw_and_save_comparison(xc, T_proj,T_vtk,H_proj,H_vtk):
    # plot the data
    file_prefix = 'verification_pure_water/'
    figure_data = {
        'T': (file_prefix + 'figure_4e_temperature_at_1500_years.png', 'T - Fig. 4E P. WEIS (2014)', 'T - VTKsample + GEOMAR'),
        'H': (file_prefix + 'figure_4e_enthalpy_at_1500_years.png', 'H - Fig. 4E P. WEIS (2014)', 'H - VTKsample + GEOMAR'),
    }
    fields_data = {
        'T': (T_proj,T_vtk),
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

def __bisection(flasher, p, s, tol=1e-8, max_iter=1000):
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

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        s_v = 1.0 - s_proj[i]
        if np.isclose(s_v, 0.0) or np.isclose(s_v, 1.0):
            PT = flasher.flash(P=pair[0]*1.0e6, T=pair[1])
            h_data.append(PT.H_mass()*1.0e-6)
        else:
            vf = __bisection(flasher, pair[0], s_v)
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

    h_data = []
    for i, pair in enumerate(zip(p_proj, t_proj)):
        s_v = 1.0 - s_proj[i]
        if np.isclose(s_v, 0.0) or np.isclose(s_v, 1.0):
            PT = flasher.flash(P=pair[0]*1.0e6, T=pair[1])
            h_data.append(PT.H_mass()*1.0e-6)
        else:
            vf = __bisection(flasher, pair[0], s_v)
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