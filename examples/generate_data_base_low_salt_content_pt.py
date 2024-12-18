import time

import numpy as np
import pyvista as pv
from globals import source_name
from globals import requested_fields
from globals import fields_to_smooth
from scipy.ndimage import gaussian_filter
from generator.source_handler import run_command

# Name association for vtk files
short_hand_name_map = {
    "saltwatereos-master": "original",
    "saltwatereos-enhanced_computations": "modified_low_salt_content",
}

# The parametric space (z,T,P) in the format: var_min / var_increment / var_max
parametric_space = {
    0: "0.00001/0.00001/0.0000201/280.0/50.0/1000.0/6/2.0/600",
    1: "0.00001/0.00001/0.0000201/280.0/25.0/1000.0/6/1.0/600",
    2: "0.00001/0.00001/0.0000201/280.0/12.5/1000.0/6/0.5/600",
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
        xhp_space = pv.read(untouched_vtk_file)
        te = time.time()
        print("Loading data: Elapsed time: ", te - tb)

        tb = time.time()
        for field in requested_fields:
            if field not in xhp_space.array_names:
                raise ValueError("The vtk file does not contain the field: ", field)

        gradients = {}
        for field in requested_fields:
            if field in ['Xv', 'Xl',  'S_h']:
                xhp_space.point_data[field] *= 0.0
            grad_field = xhp_space.compute_derivative(field, gradient=True)
            gradients[field] = grad_field["gradient"]

        for field in requested_fields:
            xhp_space.point_data.set_vectors(gradients[field], "grad_" + field)

        vtk_file_with_gradients = (
            "XTP_l" + str(level) + "_" + short_hand_name_map[source_name] + ".vtk"
        )
        xhp_space.save(vtk_file_with_gradients, binary=True)
        te = time.time()
        print("Computing gradients on requested fields: Elapsed time: ", te - tb)

    print("Process complete.")
