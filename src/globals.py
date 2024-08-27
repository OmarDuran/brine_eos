# main directive
original_source_q: bool = False


# source code
if original_source_q:
    # original source
    source_url = "https://github.com/zguoch/saltwatereos/archive/refs/heads/master.zip"
    source_name = "saltwatereos-master"
else:
    # modified source
    source_url = "https://github.com/OmarDuran/saltwatereos/archive/refs/heads/enhanced_computations.zip"
    source_name = "saltwatereos-enhanced_computations"


# requested fields
if original_source_q:
    # original source
    requested_fields = [
        "H",
        "H_l",
        "H_v",
        "H_h",
        "Rho",
        "Rho_l",
        "Rho_v",
        "Rho_h",
        "Temperature",
        "Xl",
        "Xv",
        "mu_l",
        "mu_v",
    ]
else:
    # modified source
    requested_fields = [
        "H",
        "H_l",
        "H_v",
        "H_h",
        "Rho",
        "Rho_l",
        "Rho_v",
        "Rho_h",
        "S_l",
        "S_v",
        "S_h",
        "Temperature",
        "Xl",
        "Xv",
        "mu_l",
        "mu_v",
    ]

# fields to smooth
if original_source_q:
    # original source
    fields_to_smooth = [
        "H_l",
        "H_v",
        "H_h",
        "Rho",
        "Rho_l",
        "Rho_v",
        "Rho_h",
        "Temperature",
        "Xl",
        "Xv",
        "mu_l",
        "mu_v",
    ]
else:
    # modified source
    fields_to_smooth = [
        "H_l",
        "H_v",
        "H_h",
        "Rho",
        "Rho_l",
        "Rho_v",
        "Rho_h",
        "S_l",
        "S_v",
        "S_h",
        "Temperature",
        "Xl",
        "Xv",
        "mu_l",
        "mu_v",
    ]