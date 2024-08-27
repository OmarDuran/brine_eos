import time
import pyvista as pv
from globals import source_name
from globals import requested_fields
from generator.source_handler import run_command

# Name association for vtk files
short_hand_name_map = {
    "saltwatereos-master": "original",
    "saltwatereos-enhanced_computations": "modified",
}

# The parametric space (z,T,P) in the format: var_min / var_increment / var_max
parametric_space = {
    0: "0.00001/0.1/0.511/1.0/50.0/956.84/11/50.0/620",
    1: "0.00001/0.05/0.515/1.0/25.0/956.84/11/25.0/620",
    2: "0.00001/0.025/0.5125/1.0/12.5/956.84/11/12.5/620",
}

try:
    prefix = "./cli_" + source_name + "/"
    for level, par_space in parametric_space.items():
        suffix = "XTP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        run_command(prefix + "swEOS -D 3 -V XTP -R " + par_space + " -O " + suffix)
except Exception as e:
    print(f"Error: {e} - An unexpected error occurred.")
finally:

    for level, par_space in parametric_space.items():
        untouched_vtk_file = (
            "XTP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        )
        print(" Computing gradients for file: ", untouched_vtk_file)
        # include gradients
        tb = time.time()
        xtp_space = pv.read(untouched_vtk_file)
        te = time.time()
        print("Loading data: Elapsed time: ", te - tb)

        tb = time.time()
        for field in requested_fields:
            if field not in xtp_space.array_names:
                raise ValueError("The vtk file does not contain the field: ", field)

        gradients = {}
        for field in requested_fields:
            grad_field = xtp_space.compute_derivative(field, gradient=True)
            gradients[field] = grad_field["gradient"]

        for field in requested_fields:
            xtp_space.point_data.set_vectors(gradients[field], "grad_" + field)

        vtk_file_with_gradients = (
            "XTP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        )
        xtp_space.save(vtk_file_with_gradients, binary=True)
        te = time.time()
        print("Computing gradients on requested fields: Elapsed time: ", te - tb)

    print("Process complete.")
