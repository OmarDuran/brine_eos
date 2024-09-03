import time
from water_utils import compute_flasher_data
from water_utils import generate_cartesian_grid, generate_simplex_grid_pure_water
from water_utils import assign_point_data_with_flasher_data
from water_utils import compute_gradients_and_save_vtk

H_range = (1.0e-5, 3.5)
P_range = (1.0e-3, 20.0)

for level in range(1):
    tb = time.time()
    print('Generating data base for level: ', level)
    # Generate parametric mesh
    # mesh = generate_cartesian_grid(level, H_range, P_range)
    mesh = generate_simplex_grid_pure_water(level, H_range, P_range)
    print('Number of mesh points: ', len(mesh.points))

    # Compute phase equilibrium at each mesh point
    flasher_data = compute_flasher_data(mesh.points)

    # # Assign flash fields as point data
    assign_point_data_with_flasher_data(mesh, flasher_data)

    # Compute gradients and save the data base as vtk
    compute_gradients_and_save_vtk(level, mesh)
    te = time.time()
    print("Process complete: Elapsed time: ", te - tb)
    print('')

