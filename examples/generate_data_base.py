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

# The parametric space (z,H,P) in the format: var_min / var_increment / var_max
parametric_space = {
    0: "0.0001/0.1/0.3/100/100/3000/11/100/500",
    1: "0.0001/0.05/0.3/100/50/3000/11/50/500",
    2: "0.0001/0.025/0.3/100/25/3000/11/25/500",
}

try:
    prefix = "./cli_" + source_name + "/"
    for level, par_space in parametric_space.items():
        suffix = "XHP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        run_command(prefix + "swEOS -D 3 -V XHP -R " + par_space + " -O " + suffix)
except Exception as e:
    print(f"Error: {e} - An unexpected error occurred.")
finally:

    for level, par_space in parametric_space.items():
        untouched_vtk_file = (
            "XHP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        )
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
            grad_field = xhp_space.compute_derivative(field, gradient=True)
            gradients[field] = grad_field["gradient"]

        for field in requested_fields:
            xhp_space.point_data.set_vectors(gradients[field], "grad_" + field)

        vtk_file_with_gradients = (
            "XHP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        )
        xhp_space.save(vtk_file_with_gradients, binary=True)
        te = time.time()
        print("Computing gradients on requested fields: Elapsed time: ", te - tb)

    print("Process complete.")
