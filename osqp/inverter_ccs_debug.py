import os
import pickle
import plecs

import osqp

import inverter_ccs

import scipy.signal
import scipy.sparse

import pyctl as ctl
import numpy as np
import matplotlib.pyplot as plt

from distutils.dir_util import copy_tree

import opt.osqp_gen_st

plt.ion()

# --- Input ---
export_dmpc_path = ''

ii_lim = 15.0

scaling = 200.0

t_sim = 41e-3

ref_ini = [10, 0]
ref_final = [-5, -5]
t_ref_change = 20e-3

# --- Sets up the CCS problem ---
k = 0.875
r_w = 2.e-5
n = 4
n_r = 2

file_ccs = 'model-ccs'
file_path_ccs = os.path.abspath(os.getcwd())

sim_params_ccs = plecs.Params().params
sim_params_ccs['ts'] = 1 / 7.5e3
sim_params_ccs['t_sim'] = t_sim
sim_params_ccs['i_ref_init'] = ref_ini
sim_params_ccs['i_ref_final'] = ref_final
sim_params_ccs['t_ref_step'] = t_ref_change - 1/7.5e3 - 1/80e3

ctl_params_ccs = {}
ctl_params_ccs['scaling'] = scaling

sim_params_ccs['V_dc'] = 400

V_dc = sim_params_ccs['V_dc']
u_lim = V_dc / np.sqrt(3)
u_d_lim = k * u_lim
u_q_lim = np.sqrt(1 - k ** 2) * u_lim

inv_ccs = inverter_ccs.Inverter_CCS(sim_params_ccs, ctl_params=ctl_params_ccs)
inv_ccs.set_rw(rw_1=r_w, rw_2=r_w)
inv_ccs.set_horizons(n_p=n, n_c=n, n_r=n_r)
inv_ccs.set_x_lim(ii_d=ii_lim, ii_q=ii_lim)
inv_ccs.set_u_lim(u_d=u_d_lim, u_q=u_q_lim)
#inv_ccs.set_c_export_file(export_dmpc_path)


# --- Matrices for OSQP ---
n_u = 2
n_x_ct = 2
bounds_size = n_u * n_r + n_x_ct * n_r
lin_cost_size = n_u * n

sys = inv_ccs.update_c_controller(retsys=True)
P = sys.constr_model.E_j
P = scipy.sparse.csc_matrix(P)

A = np.zeros([bounds_size, lin_cost_size])
A[:(n_u * n_r), :] = sys.constr_model.M[(n_u * n_r):2*(n_u * n_r), :]
A[(n_u * n_r):, :] = sys.constr_model.M[(2*(n_u * n_r) + n_x_ct * n_r):, :]
A = scipy.sparse.csc_matrix(A)

l = np.zeros(bounds_size)
l[:n_u * n_r] = np.tile([-u_d_lim, -u_q_lim], n_r)
l[n_u * n_r:] = np.tile([-ii_lim, -ii_lim], n_r)
l = l / scaling

u = np.zeros(bounds_size)
u[:n_u * n_r] = np.tile([u_d_lim, u_q_lim], n_r)
u[n_u * n_r:] = np.tile([ii_lim, ii_lim], n_r)
u = u / scaling

q = -sys.constr_model.Phi.T @ sys.constr_model.R_s_bar @ np.array(ref_ini).reshape(-1, 1) / scaling

prob = osqp.OSQP()
#prob.setup(P, q, A, l, u, scaled_termination=True, check_termination=40)
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=1)
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=0, max_iter=150)
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=0, max_iter=25)

#prob.setup(P, q, A, l, u, scaled_termination=True, check_termination=15, max_iter=100, warm_start=True) # 4 x 1
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=15, max_iter=100, warm_start=True) # 4 x 1
prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=30, max_iter=60, warm_start=False) # 4 x 1

#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=25, max_iter=100) # 4 x 1
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=0, max_iter=25) # 4 x 1
#prob.setup(P, q, A, l, u, scaled_termination=False, check_termination=0, max_iter=50) # 4x4

res = prob.solve()
#print(res.x)

prob.codegen('code', parameters='vectors', force_rewrite=True, FLOAT=True, LONG=False, compile_python_ext=False)
copy_tree('code/include', "C:/Users/mguerreiro/Documents/Projects/C/cdmpc/osqp")
copy_tree('code/src/osqp', "C:/Users/mguerreiro/Documents/Projects/C/cdmpc/osqp")

# --- Modified code gen ---
path_st = r'C:\Users\mguerreiro\Documents\Projects\C\cdmpc\osqp'
qdldl_fname = 'qdldl_st'
auxil_fname = 'auxil_st'

L_p_osqp = prob._model._get_workspace()['linsys_solver']['L']['p']
L_i_osqp = prob._model._get_workspace()['linsys_solver']['L']['i']
P_osqp = prob._model._get_workspace()['linsys_solver']['P']
alpha_osqp = prob._model._get_workspace()['settings']['alpha']
n_osqp = prob._model._get_workspace()['data']['n']
m_osqp = prob._model._get_workspace()['data']['m']

Dinv_osqp = prob._model._get_workspace()['linsys_solver']['Dinv']
L_x_osqp = prob._model._get_workspace()['linsys_solver']['L']['x']

#opt.osqp_gen_st.gen_qdldl_solve(qdldl_fname, L_p_osqp, L_i_osqp, P_osqp, path=path_st)

opt.osqp_gen_st.gen_qdldl_solve_v2(qdldl_fname, L_p_osqp, L_i_osqp, P_osqp, Dinv_osqp, L_x_osqp, path=path_st)
opt.osqp_gen_st.gen_auxil(auxil_fname, n_osqp, m_osqp, alpha_osqp, path=path_st)

# --- Runs PLECS sim ---
t_ccs, data_ccs = inv_ccs.run_plecs_sim(file_ccs, file_path_ccs, close=False)

ii_ccs, ig_ccs, vc_ccs, vg_ccs = inv_ccs.decode_abc(data_ccs)
ii_dq_ccs, ig_dq_ccs, vc_dq_ccs, vg_dq_ccs = inv_ccs.decode_dq0(data_ccs)

# --- Plots ---

##plt.figure()
##plt.plot(t_ccs, ig_ccs)
##plt.grid()
##
##plt.figure()
##plt.plot(t_ccs, ig_dq_ccs)
##plt.grid()

# --- Saves the sim data ---
#sim_data = [inv_ccs, inv_fcs]

##with open(data_file, 'wb') as f:
##    pickle.dump(sim_data, f)
