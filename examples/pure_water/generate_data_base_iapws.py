import time
from water_utils import compute_flasher_data
from water_utils import generate_cartesian_grid, generate_simplex_grid
from water_utils import assign_point_data_with_flasher_data
from water_utils import compute_gradients_and_save_vtk

simplex_q = True
H_range = (1.0e-5, 4.0)
P_range = (1.0e-3, 21.0)

for level in range(3):
    tb = time.time()
    print('Generating data base for level: ', level)
    # Generate parametric mesh
    if simplex_q:
        mesh = generate_simplex_grid(level, H_range, P_range)
    else:
        mesh = generate_cartesian_grid(level, H_range, P_range)

    print('Number of mesh points: ', len(mesh.points))

    # Compute phase equilibrium at each mesh point
    flasher_data = compute_flasher_data(mesh.points)

    # Assign flash fields as point data
    assign_point_data_with_flasher_data(mesh, flasher_data)

    # Compute gradients and save the data base as vtk
    compute_gradients_and_save_vtk(level, mesh, simplex_q)
    te = time.time()
    print("Process complete: Elapsed time: ", te - tb)
    print('')

