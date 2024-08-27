import time
import numpy as np
import pyvista as pv
import meshio
from globals import requested_fields
from globals import fields_to_smooth
from scipy.ndimage import gaussian_filter
from thermo import FlashPureVLS, IAPWS95Liquid, IAPWS95Gas, iapws_constants, iapws_correlations

# triple point of water
T_ref = 273.16
P_ref = 611.657
MW_H2O = iapws_constants.MWs[0] # [Kg/mol]
liquid = IAPWS95Liquid(T=T_ref, P=P_ref, zs=[1])
gas = IAPWS95Gas(T=T_ref, P=P_ref, zs=[1])
flasher = FlashPureVLS(iapws_constants, iapws_correlations, gas, [liquid], [])


ref_data={
    0: (75,75),
    1: (150,150),
    2: (300,300),
}

level = 1
H_vals = np.linspace(10, 3000.0, ref_data[level][0]) # [KJ/Kg]
P_vals = np.linspace(1.0, 50.0, ref_data[level][1]) # [MPa]

x = np.array([0.0,1.0])        # 1D array of shape (3,)
y = H_vals          # 1D array of shape (2,)
z = P_vals     # 1D array of shape (4,)

# Create the meshgrid
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

H_scale = 1.0 #1.0e-6

def S_l(rho, rho_l, rho_v):
    s_v = (rho - rho_v) / (rho_l - rho_v)
    return s_v
# Density plots
n_data = 13
flasher_data = np.empty((0, n_data), float)
for z_val, H_val, P_val in zip(X.ravel(), Y.ravel(), Z.ravel()):
    PH = flasher.flash(P=P_val*1e5, H=H_val*MW_H2O)
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
            flasher_data = np.append(flasher_data,np.array([data]),axis=0)
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
            flasher_data = np.append(flasher_data,np.array([data]),axis=0)
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
            # S_l(PH.V(), PH.phases[1].V(), PH.phases[0].V()),
            # 1.0 - S_l(PH.V(), PH.phases[1].V(), PH.phases[0].V()),
            PH.T,
            0.0,
            0.0,
            PH.phases[1].kinematic_viscosity(),
            PH.phases[0].kinematic_viscosity(),
        ]
        flasher_data = np.append(flasher_data, np.array([data]), axis=0)


# Step 2: Flatten the mesh grids and create points
points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

# Optionally create connectivity (if you need cells, e.g., for tetrahedrons, hexes, etc.)
# Here we'll just create points, but you can define cells as needed.
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

# Step 4: Assign point data (e.g., scalar field or vector field)
# Example: Assign a scalar value equal to the sum of the coordinates
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

# Save the mesh to a file
file_name = 'XHP_l'+str(level)+'_iapws_modified.vtk'

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