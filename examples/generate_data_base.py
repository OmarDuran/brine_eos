from globals import source_name
from generator.source_handler import run_command

# name association for vtk files
short_hand_name_map = {
    'saltwatereos-master': 'original',
    'saltwatereos-enhanced_computations': 'modified',
}

# The parametric space (z,H,P) in the format: var_min / var_increment / var_max
parametric_space = {
    0: '0.0001/0.1/0.3/100/100/3000/11/100/500',
    1: '0.0001/0.05/0.3/100/50/3000/11/50/500',
    # 2: '0.0001/0.025/0.3/100/25/3000/11/25/500',
}

try:
    prefix = './cli_' + source_name + "/"
    for level, par_space in parametric_space.items():
        suffix = 'XHP_l' + str(level) + '_'+ short_hand_name_map[source_name] + '.vtk'
        run_command(prefix + "swEOS -D 3 -V XHP -R " + par_space + " -O " + suffix)
except Exception as e:
    print(f"Error: {e} - An unexpected error occurred.")
finally:
    print("Process complete.")

