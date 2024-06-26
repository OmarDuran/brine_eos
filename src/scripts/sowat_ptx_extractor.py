import numpy as np
import os
from pathlib import Path
from generator.source_handler import run_command

__phase_state_map = {
    "single": "F",
    "L+H": "L+H",
    "V+L": "V+L",
    "V+H": "V+H",
}


def __retrieve_line_idx_with_pattern(
    pattern: str ,
    origin_folder: str = "apple_script/",
    origin_file_name: str = "raw_data_sowat_ptx.txt",
):
    awk_command = (
        "awk -f "
        + origin_folder
        + '/scan_line_idxs.awk ' + pattern + " "
        + origin_folder
        + origin_file_name
        + " > lines_indices.txt"
    )
    run_command(awk_command)

    idx_file_name = "lines_indices.txt"
    file = open(idx_file_name)
    line_idxs = [int(idx) - 1 for idx in file.readlines()]
    popped_item = line_idxs.pop()
    print("Deleted index: ", popped_item)
    file.close()

    clean_command = "rm -r lines_indices.txt"
    run_command(clean_command)
    return line_idxs


def __extract_PTZ_states(line_idx, content):
    T_line = content[line_idx + 0].split(" ")
    P_line = content[line_idx + 1].split(" ")
    Z_line = content[line_idx + 2].split(" ")
    return [float(Q_line[-1]) for Q_line in [P_line, T_line, Z_line]]


def __extract_phase_state(line_idx, content):
    state_line = content[line_idx + 4].split("Fluid is in the")
    state_str = "".join(
        [
            chunk
            for chunk in state_line[1].split(" ")
            if chunk not in ["", "phase", "state", "*\n"]
        ]
    )
    return __phase_state_map[state_str]


def __extract_properties(states, phase_state, line_idx, content):

    p, t, z = states
    rho_f = lambda rho_vec, S_vec: np.dot(rho_vec, S_vec)
    s_a = lambda rho_a, rho_b, Xa, Xb, z: ((Xb - z) * rho_b) / (
        -(Xa * rho_a) + z * (rho_a - rho_b) + Xb * rho_b
    )

    # extracting properties by case
    if phase_state == "F":  # essentially single phase
        rho_l_line = content[line_idx + 7]
        Xl_line = content[line_idx + 10]
        rho_l = float(rho_l_line.split(" ")[-2])
        Xl = float(Xl_line.split(" ")[-2])
        Xv = 0.0

        S_l = 1.0
        S_v = 0.0
        S_h = 0.0

        rho = rho_l  # same as rho_f(np.array([rho_l]), np.array([S_l]))
        rho_v = rho  # extended density
        rho_h = rho  # extended density
        props = [rho, rho_l, rho_v, rho_h, Xl, Xv, S_l, S_v, S_h]
        return props
    elif phase_state == "L+H":  # essentially two phase

        rho_l_line = content[line_idx + 7]
        Xl_line = content[line_idx + 10]
        rho_l = float(rho_l_line.split(" ")[-2])
        Xl = float(Xl_line.split(" ")[-2])

        rho_h_line = content[line_idx + 13]
        Xl_line = content[line_idx + 16]
        rho_h = float(rho_h_line.split(" ")[-2])

        Xh = float(Xl_line.split(" ")[-2])
        Xv = 0.0

        S_l = s_a(rho_l, rho_h, Xl, Xh, z)
        S_v = 0.0
        S_h = 1.0 - s_a(rho_l, rho_h, Xl, Xh, z)

        rho = rho_f(np.array([rho_l, rho_h]), np.array([S_l, S_h]))
        rho_v = rho  # extended density

        props = [rho, rho_l, rho_v, rho_h, Xl, Xv, S_l, S_v, S_h]
        return props
    elif phase_state == "V+L":  # essentially two phase

        rho_v_line = content[line_idx + 7]
        Xv_line = content[line_idx + 10]
        rho_v = float(rho_v_line.split(" ")[-2])
        Xv = float(Xv_line.split(" ")[-2])

        rho_l_line = content[line_idx + 13]
        Xl_line = content[line_idx + 16]
        rho_l = float(rho_l_line.split(" ")[-2])
        Xl = float(Xl_line.split(" ")[-2])

        S_l = 1.0 - s_a(rho_v, rho_l, Xv, Xl, z)
        S_v = s_a(rho_v, rho_l, Xv, Xl, z)
        S_h = 0.0

        rho = rho_f(np.array([rho_l, rho_v]), np.array([S_l, S_v]))
        rho_h = rho  # extended density

        props = [rho, rho_l, rho_v, rho_h, Xl, Xv, S_l, S_v, S_h]
        return props
    elif phase_state == "V+H":  # essentially two phase
        rho_v_line = content[line_idx + 7]
        Xv_line = content[line_idx + 10]
        rho_v = float(rho_v_line.split(" ")[-2])
        Xv = float(Xv_line.split(" ")[-2])

        rho_h_line = content[line_idx + 13]
        Xh_line = content[line_idx + 16]
        rho_h = float(rho_h_line.split(" ")[-2])
        Xh = float(Xh_line.split(" ")[-2])
        Xl = 0.0

        S_l = 0.0
        S_v = s_a(rho_v, rho_h, Xv, Xh, z)
        S_h = 1.0 - s_a(rho_v, rho_h, Xv, Xh, z)

        rho = rho_f(np.array([rho_v, rho_h]), np.array([S_v, S_h]))
        rho_l = rho  # extended density
        props = [rho, rho_l, rho_v, rho_h, Xl, Xv, S_l, S_v, S_h]
        return props
    else:
        raise ValueError("Phase state not implemented: ", phase_state)


def harvest_data(
    destination_file_name: str = "extracted_data_sowat_ptx.txt",
):

    # search for the persistent raw file (raw_data_sowat_ptx.txt)
    current_script_directory = Path(__file__).parent
    origin_folder = str(current_script_directory / "apple_script")
    origin_file_name: str = "/raw_data_sowat_ptx.txt"

    pattern = '"Please enter T in"'
    line_indices = __retrieve_line_idx_with_pattern(pattern, origin_folder, origin_file_name)

    extracted_data = np.empty((0, 12), float)
    with open(origin_folder + origin_file_name, "r") as file:
        # Read the file line by line
        content = file.readlines()

        for line_idx in line_indices:
            print("line_idx: ", line_idx)
            if line_idx == 46788:
                continue
            states = __extract_PTZ_states(line_idx, content)
            phase_state = __extract_phase_state(line_idx, content)
            properties = __extract_properties(states, phase_state, line_idx, content)
            chunk = states + properties
            extracted_data = np.append(extracted_data, np.array([chunk]), axis=0)


    normal_header = "P [bar], T [C],  z_NaCl ,  rho [Kg/m3], rho_l [Kg/m3], rho_v [Kg/m3], rho_h, Xl, Xv, S_l, S_v, S_h"
    np.savetxt(
        destination_file_name,
        extracted_data,
        delimiter=",",
        fmt="%1.6f",
        header=normal_header,
    )
