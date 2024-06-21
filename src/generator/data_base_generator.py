import os


def range_to_str(q_data):
    q_min, q_incr, q_max = q_data
    str_range = str(q_min) + "/" + str(q_incr) + "/" + str(q_max)
    return str_range

def compose_command(PHX_q, coarsening_level, p_data, e_data, x_data):
    prefix = "./swEOS -D 3 -V "
    if PHX_q:
        prefix += "PHX -R "
    else:
        prefix += "PTX -R "

    suffix = " -O PTX_c" + str(coarsening_level) + ".vtk"
    if PHX_q:
        suffix = " -O PHX_c" + str(coarsening_level) + ".vtk"

    str_range =  range_to_str(p_data)+"/"+range_to_str(e_data)+"/"+range_to_str(x_data)
    command = prefix + str_range + suffix
    return command

PHX_q = True
coarsening_level = 0

if PHX_q:
    e_min, e_incr, e_max = [2580.66,40,4640.83]
else:
    e_min, e_incr, e_max = [0.5,5,1000]

p_min, p_incr, p_max = [5,5,1000]
x_min, x_incr, x_max = [0.0,0.005,0.5]

for i, mult in enumerate(range(coarsening_level+1)):
    p_incr *= (i+1)
    e_incr *= (i+1)
    x_incr *= (i+1)

command = compose_command(PHX_q, coarsening_level,(p_min, p_incr, p_max),(e_min, e_incr, e_max),(x_min, x_incr, x_max))
print("Command: ", command)


