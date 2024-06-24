import time
import numpy as np
from sampler.vtk_sampler import VTKSampler
import matplotlib.pyplot as plt

folder_name = "Brine_isotherms_LV"
import os

if not os.path.exists(folder_name):
    os.makedirs(folder_name)


def Find_Hm(sampler_obj, T_spec, z, p):

    max_iter = 500
    tol = 1.0e-5

    # Residual function for bisection method
    def residual_T(T_spec, T_val):
        res = np.ones_like(T_val) - T_val / T_spec
        return res

    # fields initialization
    H_a = 100.0 * np.ones_like(p) * (1.0e3)
    H_b = 3000.0 * np.ones_like(p) * (1.0e3)
    H_m = 0.5 * (H_b + H_a)

    for it in range(max_iter):
        x_a = np.array((z, H_a, p)).T
        x_b = np.array((z, H_b, p)).T
        x_m = np.array((z, H_m, p)).T

        sampler_obj.sample_at(x_a)
        Ta = sampler_obj.sampled_could.point_data["Temperature"]

        sampler_obj.sample_at(x_b)
        Tb = sampler_obj.sampled_could.point_data["Temperature"]

        sampler_obj.sample_at(x_m)
        Tm = sampler_obj.sampled_could.point_data["Temperature"]

        res_a = residual_T(T_spec, Ta)
        res_b = residual_T(T_spec, Tb)
        res_m = residual_T(T_spec, Tm)

        check_points = res_a * res_b < 0
        valid_idx = np.where(check_points)[0]

        if np.all(
            np.logical_or(
                np.abs(res_m[valid_idx]) < tol,
                np.abs(H_b[valid_idx] - H_a[valid_idx]) < tol,
            )
        ):
            print("Find_Hm:: Method converge.")
            print(
                "Find_Hm:: Percentage of recovered points: ",
                100.0 * float((valid_idx.shape[0] / H_m.shape[0])),
            )
            print("Find_Hm:: Number of global iterations: ", it)
            print(
                "Find_Hm:: Overall in Temperature [K] l2_error: ",
                np.linalg.norm(Tm - T_spec),
            )
            return H_m, valid_idx

        # new states
        idx_n = np.where(res_a * res_m < 0)
        idx_p = np.where(res_a * res_m >= 0)
        H_b[idx_n] = H_m[idx_n]
        H_a[idx_p] = H_m[idx_p]
        H_m = 0.5 * (H_b + H_a)


to_K = lambda T: T + 273.15
T_spec_to_files = {
    325: "sowat_input/TD_LV_T_325.txt",
    370: "sowat_input/TD_LV_T_370.txt",
    400: "sowat_input/TD_LV_T_400.txt",
}

for level in [0, 1, 2]:
    # Create VTKSampler object from the generated vtk file
    vtk_file_name = "XHP_l" + str(level) + "_modified.vtk"
    taylor_extended = True
    sampler_obj = VTKSampler(vtk_file_name, taylor_extended)
    sampler_obj.conversion_factors = (0.01, 1.0e-3, 1.0)  # (z[%], h [kJ/kg], p [bar])

    T_spec_to_zHp = {}
    for T_spec, file_name in T_spec_to_files.items():

        file_name = T_spec_to_files[T_spec]
        # Load the data from sowatflinc_plotdata.exe
        # see https://mineralsystems.ethz.ch/software/sowat.html
        isotherm_z_p = np.loadtxt(file_name)
        z, p = isotherm_z_p.T
        tb = time.time()
        H_m, idx = Find_Hm(sampler_obj, to_K(T_spec), z, p)
        te = time.time()
        print("Overall time for computing Enthalpy: ", te - tb)
        T_spec_to_zHp[T_spec] = np.vstack((z, H_m, p)).T

    # Create a figure and axes
    fig, ax = plt.subplots()

    # Plot the lines with labels
    ax.plot(
        T_spec_to_zHp[325][:, 1],
        T_spec_to_zHp[325][:, 2],
        label="T = 325 C",
        color="r",
        linestyle="-",
        linewidth=2,
    )
    ax.plot(
        T_spec_to_zHp[370][:, 1],
        T_spec_to_zHp[370][:, 2],
        label="T = 370 C",
        color="g",
        linestyle="-",
        linewidth=2,
    )
    ax.plot(
        T_spec_to_zHp[400][:, 1],
        T_spec_to_zHp[400][:, 2],
        label="T = 400 C",
        color="b",
        linestyle="-",
        linewidth=2,
    )

    # Add axis labels and a title
    ax.set_xlabel("Mixture enthalpy [J/Kg]")
    ax.set_ylabel("Pressure [bar]")
    ax.set_title("Isotherms for refinement level: " + str(level))
    ax.legend()
    fig_temp = folder_name + "/" + "isotherms_from_sowatflinc_l" + str(level) + ".png"
    plt.savefig(fig_temp)
    # plt.show()
