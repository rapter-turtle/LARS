#!/usr/bin/env python3
"""
Generate acados C code for the HERON USV MPC model.
This script:
  - builds the OCP
  - exports solver C code into ./c_generated_code
"""

import os
import numpy as np

# 🔴 핵심 추가
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)


from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
from casadi import SX, vertcat, sin, cos, tan, sqrt, exp, tanh
import scipy.linalg


# ===========================================================
# 1) MODEL DEFINITION
# ===========================================================
def export_heron_model():
    model_name = "heron"

    # States
    xn   = SX.sym('xn')
    yn   = SX.sym('yn')
    psi  = SX.sym('psi')
    u    = SX.sym('u')
    v    = SX.sym('v')
    r    = SX.sym('r')
    delta = SX.sym('delta')
    F     = SX.sym('F')

    states = vertcat(xn, yn, psi, u, v, r, delta, F)

    # Controls
    delta_d = SX.sym('delta_d')
    F_d     = SX.sym('F_d')
    inputs  = vertcat(delta_d, F_d)


    # double obs[HERON_NP] = {
    #     CD_x, CD_y, CD_psi,
    #     MS_x, MS_y, MS_psi,
    #     dob_state_.param_filtered[0],  dob_state_.param_filtered[1],  dob_state_.param_filtered[2],
    #     stream_speed, switch_1
    # };
    # Disturbances + obstacle parameters
    CD_x = SX.sym('CD_x')
    CD_y = SX.sym('CD_y')
    CD_psi = SX.sym('or1')
    MS_x = SX.sym('ox2')
    MS_y = SX.sym('oy2')
    MS_psi = SX.sym('or2')
    du  = SX.sym('du')
    dv  = SX.sym('dv')
    dr  = SX.sym('dr')

    switch  = SX.sym('switch')
    stream_speed  = SX.sym('stream_speed')

    Xu  = SX.sym('Xu')
    Xuu = SX.sym('Xuu')
    Yv  = SX.sym('Yv')
    Yvv = SX.sym('Yvv')
    Yr  = SX.sym('Yr')
    Yrr  = SX.sym('Yrr')
    Yrrr  = SX.sym('Yrrr')
    Nr  = SX.sym('Nr')
    Nrr = SX.sym('Nrr')
    Nrrr = SX.sym('Nrrr')
    Nv  = SX.sym('Nv')
    u_F = SX.sym('u_F')
    v_F = SX.sym('v_F')
    r_F = SX.sym('r_F')
    
    Xu_dot = SX.sym('Xu_dot')
    Yv_dot = SX.sym('Yv_dot')
    Nr_dot = SX.sym('Nr_dot')
    M = SX.sym('M')
    I = SX.sym('I')    

    alpha_side = SX.sym('alpha_side')
    alpha_dock = SX.sym('alpha_dock')

    cbf_d = SX.sym('cbf_d')
    cbf_a = SX.sym('cbf_a')
    cbf_b = SX.sym('cbf_b')
    cbf_x0 = SX.sym('cbf_x0')
    cbf_y0 = SX.sym('cbf_y0')   

    p = vertcat(CD_x, CD_y, CD_psi,
                MS_x, MS_y, MS_psi,
                du, dv, dr,
                stream_speed, switch, 
                Xu, Xuu,
                Yv, Yvv, Yr, Yrr, Yrrr,
                Nr, Nrr, Nrrr, Nv, 
                u_F, v_F, r_F, 
                Xu_dot, Yv_dot, Nr_dot,
                M, I,
                alpha_side, alpha_dock,
                cbf_d, cbf_a, cbf_b, cbf_x0, cbf_y0)

    # xdot symbols
    xn_dot = SX.sym('xn_dot')
    yn_dot = SX.sym('yn_dot')
    psi_dot= SX.sym('psi_dot')
    u_dot  = SX.sym('u_dot')
    v_dot  = SX.sym('v_dot')
    r_dot  = SX.sym('r_dot')
    d_dot  = SX.sym('delta_dot')
    F_dot  = SX.sym('F_dot')

    states_dot = vertcat(xn_dot, yn_dot, psi_dot, u_dot, v_dot, r_dot, d_dot, F_dot)

    # Parameters for dynamics
    # Xu  = 0.0845
    # Xuu = 0.0195
    # Yv  = 0.0485
    # Yvv = 0.0988
    # Yr  = 0.151
    # Nr  = 0.6939
    # Nrr = 0.0
    # Nv  = 0.0
    # alpha1 = 0.0452
    # alpha2 = 2.0/500.0
    # alpha3 = 0.0188
    # alpha4 = 0.0193




    eps = 1e-5

    # Thrust model
    T = F
    hp = 3.0

    ur = u + stream_speed*cos(psi)
    vr = v - stream_speed*sin(psi)

    # Dynamics (expl)
    f_expl = vertcat(
        u*cos(psi) - v*sin(psi),# - hp*r*sin(psi),
        u*sin(psi) + v*cos(psi), ##+ hp*r*cos(psi),
        r,
        (Xu*ur + Xuu*sqrt(u*u+eps)*u + u_F*T*cos(0.01*delta))/(M+Xu_dot) - du,
        (Yv*vr + Yr*r + Yvv*sqrt(vr*vr+eps)*vr + Yrr*sqrt(r*r+eps)*r + Yrrr*r*r*r + v_F*T*sin(0.01*delta))/(M+Yv_dot) - dv,
        (Nr*r + Nv*vr + Nrr*sqrt(r*r+eps)*r + Nrrr*r*r*r - r_F*T*sin(0.01*delta))/(I+Nr_dot) - dr,
        delta_d,
        F_d
    )

    f_impl = states_dot - f_expl

    # Dummy constraint h(x)
    # alpha_side = 0.1
    # alpha_dock = 0.5

    # cbf_d = 4
    # cbf_a = 3.0
    # cbf_b = 0.5
    # cbf_x0 = -10.0
    # cbf_y0 = 3.5


    x_prime = MS_x + cbf_d*cos(MS_psi + 0.5*3.141592)
    y_prime = MS_y + cbf_d*sin(MS_psi + 0.5*3.141592)
    h_side = yn - tan(MS_psi)*xn - y_prime + tan(MS_psi)*x_prime 
    h_side_dot = u*sin(psi) + v*cos(psi) - tan(MS_psi)*(u*cos(psi) - v*sin(psi)- stream_speed)
    
    xh = xn + hp*cos(psi)    
    yh = yn + hp*sin(psi)
    xr = xh - CD_x    
    yr = yh - CD_y
    xh_dot = u*cos(psi) - v*sin(psi) - hp*r*sin(psi) - stream_speed
    yh_dot = u*sin(psi) + v*cos(psi) + hp*r*cos(psi)

    h_dock_left = cbf_a*tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr))) + cbf_y0 - sin(CD_psi)*xr + cos(CD_psi)*yr
    h_dock_right = cbf_a*tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr))) + cbf_y0 + sin(CD_psi)*xr - cos(CD_psi)*yr
    h_dock_left_dot = cbf_a*cbf_b*(- (cos(CD_psi)*xh_dot + sin(CD_psi)*yh_dot))*(1 - tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))*tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))) - sin(CD_psi)*xh_dot + cos(CD_psi)*yh_dot
    h_dock_right_dot = cbf_a*cbf_b*(- (cos(CD_psi)*xh_dot + sin(CD_psi)*yh_dot))*(1 - tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))*tanh(cbf_b*(cbf_x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))) + sin(CD_psi)*xh_dot - cos(CD_psi)*yh_dot


    h_expr = SX.zeros(3, 1)
    h_expr[0] = -(h_side_dot + alpha_side*h_side)
    h_expr[1] = switch*(h_dock_left_dot + alpha_dock*h_dock_left)
    h_expr[2] = switch*(h_dock_right_dot + alpha_dock*h_dock_right)
    # h_expr[1] = switch*(h_dock_left)
    # h_expr[2] = switch*(h_dock_right)


    model = AcadosModel()
    model.con_h_expr = h_expr
    model.p = p 
    model.f_impl_expr = f_impl
    model.f_expl_expr = f_expl
    model.x = states
    model.xdot = states_dot
    model.u = inputs
    model.name = model_name

    # store meta information
    model.x_labels = ['$x$ [m]', '$y$ [m]',  '$psi$ [rad]',  '$u$ [m/s]', '$v$ [m/s]', '$r$ [rad/s]', '$delta$ [N]', '$F$ [N]']
    model.u_labels = ['$n_1_d$ [N/s]', '$n_2_d$ [N/s]']
    model.t_label = '$t$ [s]'


    return model


# ===========================================================
# 2) OCP SETUP
# ===========================================================
def setup_trajectory_tracking(x0, N_horizon, Tf):

    ocp = AcadosOcp()

    model = export_heron_model()
    ocp.model = model

    nx = model.x.rows()
    nu = model.u.rows()
    ny = nx + nu
    ny_e = nx

    ocp.dims.N = N_horizon

    # set cost module
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.cost.cost_type_e = 'NONLINEAR_LS'

    # Q_mat = 2*np.diag([1e2, 1e2, 1e5, 1e1, 1e1, 1e1, 1e1, 1e-1])
    # Q_mat_term = 2*np.diag([1e5, 1e5, 1e4, 1e-2, 1e-1, 1e1, 1e1, 1e-1])
    # R_mat = 2*np.diag([1e2, 1e-2])

    Q_mat = 2*np.diag([1e3, 1e3, 1e5, 1e1, 1e3, 1e2, 1e1, 1e-1])
    Q_mat_term = 2*np.diag([1e5, 1e5, 1e4, 1e-2, 1e-1, 1e2, 1e1, 1e-1])
    R_mat = 2*np.diag([1e2, 1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat_term

    ocp.model.cost_y_expr = vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    ocp.constraints.x0 = x0

                # Xu, Xuu,
                # Yv, Yvv, Yr, Yrr, Yrrr,
                # Nr, Nrr, Nrrr, Nv, 
                # u_F, v_F, r_F)


    ocp.parameter_values = np.array([0.0, 0.0, 0.0, 
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     1.0, 1.0, 
                                     0.0, 0.0,
                                     0.0, 0.0, 0.0, 0.0, 0.0])

    num_obs = 3
    ocp.constraints.uh = 1e10 * np.ones(num_obs)
    ocp.constraints.lh = np.zeros(num_obs)


    ocp.constraints.idxsh = np.array([0,1,2])
    ocp.constraints.idxsh_e = np.array([0,1,2])
    Zh = 1e5 * np.ones(num_obs)
    zh = 1e5 * np.ones(num_obs)
    ocp.cost.zl = zh
    ocp.cost.zu = zh
    ocp.cost.Zl = Zh
    ocp.cost.Zu = Zh
    ocp.cost.zl_e = zh
    ocp.cost.zu_e = zh
    ocp.cost.Zl_e = Zh
    ocp.cost.Zu_e = Zh

    # copy for terminal
    ocp.constraints.uh_e = ocp.constraints.uh
    ocp.constraints.lh_e = ocp.constraints.lh
    ocp.model.con_h_expr_e = ocp.model.con_h_expr

    # set constraints
    steer_max = 100*60.0*3.14192/180.0
    ocp.constraints.lbu = np.array([-steer_max*0.3,-10.0])
    ocp.constraints.ubu = np.array([+steer_max*0.3,+10.0])
    ocp.constraints.idxbu = np.array([0, 1])

    ocp.constraints.lbx = np.array([-steer_max, -5.0])
    ocp.constraints.ubx = np.array([steer_max, 20])    
    ocp.constraints.idxbx = np.array([6, 7])



    # Solver options

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.sim_method_newton_iter = 150
    ocp.solver_options.nlp_solver_type = 'SQP_RTI'
    ocp.solver_options.qp_solver_cond_N = N_horizon

    
    ocp.solver_options.tf = Tf


    ocp.code_export_directory = "c_generated_code"


    return ocp


# ===========================================================
# 3) GENERATE C CODE
# ===========================================================
def main():

    print("→ Building HERON OCP model...")

    x0 = np.array([0, 0, 0, 0, 0, 0, 0, 0])
    N = 40
    Tf = 10.0

    ocp = setup_trajectory_tracking(x0, N, Tf)

    # 현재 acados_template 버전에 맞는 방식: generate=True
    print("Generating solver + exporting C code...")

    # 현재 acados_template 버전에 맞는 방식: generate=True
    solver = AcadosOcpSolver(
        ocp,
        json_file="acados_ocp_heron.json",
        generate=True,
    )

    out_dir = os.path.join(os.path.dirname(__file__), "c_generated_code")
    print(f"✔ DONE. C code exported to: {out_dir}")




if __name__ == "__main__":
    main()

