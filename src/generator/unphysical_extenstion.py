import numpy as np
import pyvista
import time

level = 0
file_name_step_1 = 'XHP_l'+ str(level)+ '.vtk'
file_name_step_2 = 'PHX_l'+str(level)+'_with_saturation.vtk'
file_name_step_3 = 'PHX_l'+str(level)+'_with_gradients.vtk'

tb = time.time()
phx_space = pyvista.read(file_name_step_1)
te = time.time()
print("Loading data: Elapsed time: ", te-tb)


# Compute and append saturation 7 possible phase states
# def define_cases(rho_chunk):
#     # rho_l, rho_v, rho_h = rho_chunk
#     nz_idx = np.nonzero(rho_chunk)
#     if len(nz_idx) == 1:
#         if nz_idx[0][0] == 0:
#             return 0 # Single phase(Liquid)
#         elif nz_idx[0][0] == 1:
#             return 1 # Single phase(Vapor)
#         else:
#             return 2 # Single phase(Halite)
#     elif len(nz_idx) == 2:
#         if nz_idx[0][0] == 0 and nz_idx[0][1] == 1:
#             return 3 # Two phase(Liquid-Vapor)
#         elif nz_idx[0][0] == 0 and nz_idx[0][1] == 2:
#             return 4 # Two phase(Liquid-Halite)
#         else:
#             return 5  # Two phase(Vapor-Halite)
#     else:
#         return 6 # Three phase(Liquid-Vapor-Halite)

def s_alpha(rho_alpha, rho_beta, rho):
    den = (rho_alpha - rho_beta)
    num = (rho - rho_beta)
    return num / den

def s_alpha_beta(rho_chunk, h_chunk):
    rho_alpha, rho_beta, rho_gamma, rho = rho_chunk
    h_alpha, h_beta, h_gamma, h = h_chunk
    den = h_alpha * rho_alpha *(rho_beta - rho_gamma) + h_gamma * (rho_alpha - rho_beta) * rho_gamma + h_beta * rho_beta * (-rho_alpha + rho_gamma)
    num_alpha = h * rho * (rho_beta - rho_gamma) +  h_gamma * (rho - rho_beta) * rho_gamma +  h_beta * rho_beta * (-rho + rho_gamma)
    num_beta = h_alpha * rho_alpha * (rho - rho_gamma) + h_gamma * (- rho + rho_alpha) * rho_gamma +  h * rho * (-rho_alpha + rho_gamma)
    s_alpha = num_alpha / den
    s_beta = num_beta / den
    s_gamma = 1.0 - s_alpha - s_beta
    return np.array([s_alpha, s_beta, s_gamma])

def compute_saturation(chunk):
    rho_chunk = chunk[np.array([0, 1, 2, 3])]
    h_chunk = chunk[np.array([4, 5, 6, 7])]
    rho_l, rho_v, rho_h, rho = rho_chunk
    nz_idx = np.nonzero(np.array([rho_l, rho_v, rho_h]))
    if len(nz_idx[0]) == 1:
        if nz_idx[0][0] == 0:
            return np.array([1.0,0.0,0.0]) # Single phase(Liquid)
        elif nz_idx[0][0] == 1:
            return np.array([0.0,1.0,0.0]) # Single phase(Vapor)
        else:
            return np.array([0.0,0.0,1.0]) # Single phase(Halite)
    elif len(nz_idx[0]) == 2:
        if nz_idx[0][0] == 0 and nz_idx[0][1] == 1:
            s_l = s_alpha(rho_l, rho_v, rho)
            return np.array([s_l,1.0-s_l,0.0]) # Two phase(Liquid-Vapor)
        elif nz_idx[0][0] == 0 and nz_idx[0][1] == 2:
            s_l = s_alpha(rho_l, rho_h, rho)
            return np.array([s_l,0.0,1.0-s_l]) # Two phase(Liquid-Halite)
        else:
            s_v = s_alpha(rho_v, rho_h, rho)
            return np.array([0.0,s_v,1.0-s_v])  # Two phase(Vapor-Halite)
    else:
        return s_alpha_beta(rho_chunk, h_chunk) # Three phase(Liquid-Vapor-Halite)

def extend_rho_and_h(chunk):
    rho_chunk = chunk[np.array([0, 1, 2, 3])]
    h_chunk = chunk[np.array([4, 5, 6, 7])]
    rho_l, rho_v, rho_h, rho = rho_chunk
    h_l, h_v, h_h, h = h_chunk
    nz_idx = np.nonzero(np.array([rho_l, rho_v, rho_h]))
    if len(nz_idx[0]) == 1:
        if nz_idx[0][0] == 0:
            return np.array([rho_l,rho,rho,h_l,h,h]) # Single phase(Liquid)
        elif nz_idx[0][0] == 1:
            return np.array([rho,rho_v,rho,h,h_v,h]) # Single phase(Vapor)
        else:
            return np.array([rho,rho,rho_h,h,h,h_h]) # Single phase(Halite)
    elif len(nz_idx[0]) == 2:
        if nz_idx[0][0] == 0 and nz_idx[0][1] == 1:
            s_l = s_alpha(rho_l, rho_v, rho)
            return np.array([rho_l,rho_v,rho,h_l,h_v,h]) # Two phase(Liquid-Vapor)
        elif nz_idx[0][0] == 0 and nz_idx[0][1] == 2:
            s_l = s_alpha(rho_l, rho_h, rho)
            return np.array([rho_l,rho,rho_h,h_l,h,h_h]) # Two phase(Liquid-Halite)
        else:
            s_v = s_alpha(rho_v, rho_h, rho)
            return np.array([rho,rho_v, rho_h,h,h_v, h_h])  # Two phase(Vapor-Halite)
    else:
        return np.array([rho_l, rho_v, rho_h,h_l, h_v, h_h]) # Three phase(Liquid-Vapor-Halite)

def extend_mu(mu_chunk):
    mu_l, mu_v = mu_chunk
    nz_idx = np.nonzero(np.array([mu_l, mu_v]))
    if len(nz_idx[0]) == 1:
        if nz_idx[0][0] == 0:
            return np.array([mu_l,mu_l]) # Single phase(Liquid)
        else:
            return np.array([mu_v,mu_v]) # Single phase(Vapor)
    else:
        return np.array([mu_l, mu_v]) # Two phase(Liquid-Vapor)

pt_data = phx_space.point_data
rho_h_data = np.vstack((pt_data['Rho_l'],pt_data['Rho_v'],pt_data['Rho_h'],pt_data['Rho'],pt_data['H_l'],pt_data['H_v'],pt_data['H_h'],pt_data['H'])).T
saturation_data = np.array(list(map(compute_saturation,rho_h_data)))
phx_space.point_data.set_array(saturation_data[:,0], 'S_l')
phx_space.point_data.set_array(saturation_data[:,1], 'S_v')
phx_space.point_data.set_array(saturation_data[:,2], 'S_h')
extended_rho_h_data = np.array(list(map(extend_rho_and_h,rho_h_data)))
phx_space.point_data.set_array(extended_rho_h_data[:,0],'Rho_l')
phx_space.point_data.set_array(extended_rho_h_data[:,1],'Rho_v')
phx_space.point_data.set_array(extended_rho_h_data[:,2],'Rho_h')
phx_space.point_data.set_array(extended_rho_h_data[:,3],'H_l')
phx_space.point_data.set_array(extended_rho_h_data[:,4],'H_v')
phx_space.point_data.set_array(extended_rho_h_data[:,5],'H_h')
phx_space.point_data.set_array(1.0/extended_rho_h_data[:,0],'nu_l')
phx_space.point_data.set_array(1.0/extended_rho_h_data[:,1],'nu_v')
phx_space.point_data.set_array(1.0/extended_rho_h_data[:,2],'nu_h')
mu_data = np.vstack((pt_data['mu_l'],pt_data['mu_v'])).T
extended_mu_data = np.array(list(map(extend_mu,mu_data)))
phx_space.point_data.set_array(extended_mu_data[:,0],'mu_l')
phx_space.point_data.set_array(extended_mu_data[:,1],'mu_v')

# convert T c to K
phx_space.point_data.set_array(pt_data['Temperature']+273.15,'Temperature')

phx_space.save(file_name_step_2, binary=True)

# compute gradients
request_fields = ['Temperature', 'Rho', 'H', 'Xl', 'Xv', 'Rho_l', 'Rho_v', 'Rho_h', 'H_l', 'H_v', 'H_h', 'mu_l', 'mu_v','S_l','S_v','S_h','nu_l', 'nu_v', 'nu_h']
tb = time.time()
phx_space = pyvista.read(file_name_step_2)
te = time.time()
print("Loading data: Elapsed time: ", te-tb)

fields = phx_space.array_names
gradients = {}
for field in fields:
    if field not in request_fields:
        continue
    grad_field = phx_space.compute_derivative(field, gradient=True)
    gradients[field] = grad_field['gradient']

for field in fields:
    if field not in request_fields:
        continue
    phx_space.point_data.set_vectors(gradients[field], 'grad_'+field)
phx_space.save(file_name_step_3, binary=True)

# for n_data in [5, 10, 100, 200]:
#     n_points = n_data*n_data*n_data
#     points = np.array([[150.0, 2000.0, 0.25] for i in range(n_points)])
#     point_cloud = pyvista.PolyData(points)
#     tb = time.time()
#     evaluated_could = point_cloud.sample(phx_space)
#     te = time.time()
#     print("n_points: ", n_points)
#     print("Sampling data: Elapsed time: ", te-tb)
