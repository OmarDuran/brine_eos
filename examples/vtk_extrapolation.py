import numpy as np
import pyvista as pv
from sampler.vtk_sampler import VTKSampler
import matplotlib.pyplot as plt
from matplotlib import cm

fig = plt.figure(figsize=plt.figaspect(1.0))
ax = fig.add_subplot(111, projection='3d')
ax.set_zlim(-0.025, +0.025)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("Function value")
# Set the view angle (elevation, azimuth)
ax.view_init(elev=50, azim=-60)

vtk_file_name = "sample_with_gradients.vtk"
folder_name = 'sampler_figures_3d'
import os
if not os.path.exists(folder_name):
    os.makedirs(folder_name)


def compose_figure_name(folder_name, z_val, suffix):
    fig_name = folder_name + '/'
    fig_name += 'z_' + str(z_val)
    fig_name += suffix
    return fig_name

def f(x, y, z):
    xs = 0.75 * x**4 - x**2 + 0.3
    ys = 0.75 * y**4 - y**2 + 0.3
    zs = 0.75 * z**4 - z**2 + 0.3
    return xs * ys * zs

def generate_vtk_from_function(f, vtk_file_name):

    # Create the 3D NumPy array of spatially referenced data again.
    nc = 30
    nc_x = nc
    nc_y = nc
    nc_z = nc
    dx = float(2.0 / nc_x)
    dy = float(2.0 / nc_y)
    dz = float(2.0 / nc_z)
    x = np.linspace(-1.0, 1.0, nc_x + 1)
    y = np.linspace(-1.0, 1.0, nc_y + 1)
    z = np.linspace(-1.0, 1.0, nc_z + 1)
    gxv, gyv, gzv = np.meshgrid(x, y, z)
    grid = pv.ImageData()
    grid.dimensions = gxv.shape

    # Edit the spatial reference
    grid.origin = (-1.0, -1.0, -1.0)  # The bottom left corner of the data set
    grid.spacing = (dx, dy, dz)  # These are the cell sizes along each axis
    grid.point_data["f"] = f(gxv, gyv, gzv).flatten(order="F")  # Flatten the array!


    fields = grid.array_names
    gradients = {}
    for field in fields:
        grad_field = grid.compute_derivative(field, gradient=True)
        gradients[field] = grad_field["gradient"]

    for field in fields:
        grid.point_data.set_vectors(gradients[field], "grad_" + field)
    grid.save(vtk_file_name, binary=True)


# step 1:
# generate a vtk file with a external function values
generate_vtk_from_function(f,vtk_file_name)

# step 2:
# create VTKSampler object from the generated vtk file
taylor_extended = False
sampled_obj = VTKSampler(vtk_file_name, taylor_extended)

# step 3:
# Compare exact and interpolated function at z_val
z_val = 0.0

# Evaluate function in cartesian manner
nc = 30
nc_x = nc
nc_y = nc
dx = float(2.0 / nc_x)
dy = float(2.0 / nc_y)
x = np.linspace(-1.0, 1.0, nc_x + 1)
y = np.linspace(0.0, 1.0, nc_y + 1)
gxv, gyv = np.meshgrid(x, y)
gzv = z_val * np.ones_like(x)
fe_v = f(gxv, gyv, gzv)


# plot for exact function
surf_1 = ax.plot_surface(gxv, gyv, fe_v, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Generate data for sampling
dxi = 0.1
xi = np.arange(-1.1, 0.0 + dxi, dxi)
eta = np.arange(-1.1, 1.1 + dxi, dxi)
pxv = xi
pyv = eta
pxv, pyv = np.meshgrid(pxv, pyv)
pxv_shape = pxv.shape
pzv = z_val * np.ones_like(pxv)

# Perform sampling
par_points = np.array((pxv.flatten(), pyv.flatten(), pzv.flatten())).T
sampled_obj.sample_at(par_points)
fv = sampled_obj.sampled_could.point_data["f"].reshape(pxv_shape)

# plot for interpolated function
norm = plt.Normalize(fe_v.min(), fv.max())
face_colors = cm.coolwarm(np.ones_like(norm(fv)))
rcount, ccount, _ = face_colors.shape
surf_2 = ax.plot_surface(pxv, pyv, fv, rcount=rcount, ccount=ccount,facecolors=face_colors, shade=False, linewidth=1.0)
surf_2.set_facecolor((0,0,0,0))

# compose figure name and save it
fig_temp = compose_figure_name(folder_name, z_val, '_smooth_function_3d.png')
plt.savefig(fig_temp)
plt.show()
