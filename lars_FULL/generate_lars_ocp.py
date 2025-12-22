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

    p = vertcat(CD_x, CD_y, CD_psi,
                MS_x, MS_y, MS_psi,
                du, dv, dr,
                stream_speed, switch)

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
    Xu = 0.0845
    Xuu = 0.0195
    Yv = 0.0485
    Yvv = 0.0988
    Yr = 0.151
    Nr = 0.6939
    Nrr = 0.0
    Nv = 0.0
    alpha1 = 0.0452
    alpha2 = 2.0/500.0
    alpha3 = 0.0188
    alpha4 = 0.0193

    eps = 1e-5

    # Thrust model
    T = F
    hp = 3.0

    # Dynamics (expl)
    f_expl = vertcat(
        -stream_speed + u*cos(psi) - v*sin(psi),# - hp*r*sin(psi),
        u*sin(psi) + v*cos(psi), ##+ hp*r*cos(psi),
        r,
        (-Xu*u - Xuu*sqrt(u*u+eps)*u + alpha1*T*cos(alpha2*delta)) - du,
        (-Yv*v - Yr*r - Yvv*sqrt(v*v+eps)*v + alpha3*T*sin(alpha2*delta)) - dv,
        (-Nr*r - Nv*v - Nrr*sqrt(r*r+eps)*r - alpha4*T*sin(alpha2*delta)) - dr,
        delta_d,
        F_d
    )

    f_impl = states_dot - f_expl

    # Dummy constraint h(x)
    alpha_side = 0.1
    alpha_dock = 0.5

    d = 4
    x_prime = MS_x + d*cos(MS_psi + 0.5*3.141592)
    y_prime = MS_y + d*sin(MS_psi + 0.5*3.141592)
    h_side = yn - tan(MS_psi)*xn - y_prime + tan(MS_psi)*x_prime 
    h_side_dot = u*sin(psi) + v*cos(psi) - tan(MS_psi)*(-stream_speed + u*cos(psi) - v*sin(psi))
    

    
    a = 3.0
    b = 0.5
    x0 = -10.0
    y0 = 3.5
    xh = xn + hp*cos(psi)    
    yh = yn + hp*sin(psi)
    xr = xh - CD_x    
    yr = yh - CD_y
    xh_dot = -stream_speed + u*cos(psi) - v*sin(psi) - hp*r*sin(psi)
    yh_dot = u*sin(psi) + v*cos(psi) + hp*r*cos(psi)

    h_dock_left = a*tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr))) + y0 - sin(CD_psi)*xr + cos(CD_psi)*yr
    h_dock_right = a*tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr))) + y0 + sin(CD_psi)*xr - cos(CD_psi)*yr
    h_dock_left_dot = a*b*(- (cos(CD_psi)*xh_dot + sin(CD_psi)*yh_dot))*(1 - tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))*tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))) - sin(CD_psi)*xh_dot + cos(CD_psi)*yh_dot
    h_dock_right_dot = a*b*(- (cos(CD_psi)*xh_dot + sin(CD_psi)*yh_dot))*(1 - tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))*tanh(b*(x0 - (cos(CD_psi)*xr + sin(CD_psi)*yr)))) + sin(CD_psi)*xh_dot - cos(CD_psi)*yh_dot


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

    Q_mat = 2*np.diag([1e3, 1e3, 1e-3, 1e-2, 0.0, 0.0, 1e-1, 1e0])
    R_mat = 2*np.diag([1e0, 1e-2])

    ocp.cost.W = scipy.linalg.block_diag(Q_mat, R_mat)
    ocp.cost.W_e = Q_mat

    ocp.model.cost_y_expr = vertcat(model.x, model.u)
    ocp.model.cost_y_expr_e = model.x
    ocp.cost.yref  = np.zeros((ny, ))
    ocp.cost.yref_e = np.zeros((ny_e, ))

    ocp.constraints.x0 = x0


    ocp.parameter_values = np.array([0.0, 0.0, 0.0, 
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0, 0.0,
                                     0.0, 0.0])

    num_obs = 3
    ocp.constraints.uh = 1e10 * np.ones(num_obs)
    ocp.constraints.lh = np.zeros(num_obs)


    ocp.constraints.idxsh = np.array([0,1,2])
    ocp.constraints.idxsh_e = np.array([0,1,2])
    Zh = 1e7 * np.ones(num_obs)
    zh = 1e7 * np.ones(num_obs)
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
    ocp.constraints.lbu = np.array([-200,-10.0])
    ocp.constraints.ubu = np.array([+200,+10.0])
    ocp.constraints.idxbu = np.array([0, 1])

    ocp.constraints.lbx = np.array([-250, 0.0])
    ocp.constraints.ubx = np.array([250, 40])    
    ocp.constraints.idxbx = np.array([6, 7])



    # Solver options

    ocp.solver_options.qp_solver = 'PARTIAL_CONDENSING_HPIPM' # FULL_CONDENSING_QPOASES
    ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
    ocp.solver_options.integrator_type = 'IRK'
    ocp.solver_options.sim_method_newton_iter = 100
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
    Tf = 20.0

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

