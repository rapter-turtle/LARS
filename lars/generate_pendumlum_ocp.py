#!/usr/bin/env python3
# generate_pendulum_ocp.py

import os
import warnings

# --- (옵션) 이 스크립트가 있는 폴더를 작업 디렉토리로 강제 세팅 ---
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# CasADi 버전 경고 숨기기
warnings.filterwarnings(
    "ignore",
    message="Full featured acados requires CasADi version >=",
    category=UserWarning,
)

from acados_template import AcadosOcp, AcadosOcpSolver, AcadosModel
from casadi import SX, vertcat, sin, cos
import numpy as np


def export_pendulum_model() -> AcadosModel:
    model_name = "pendulum"

    # states
    x     = SX.sym('x')      # cart position
    v     = SX.sym('v')      # cart velocity
    theta = SX.sym('theta')  # pendulum angle (theta = 0: upright)
    omega = SX.sym('omega')  # angular velocity

    x_sym = vertcat(x, v, theta, omega)

    # input
    F = SX.sym('F')          # horizontal force on cart
    u_sym = vertcat(F)

    # parameters
    M = 1.0    # cart mass
    m = 0.1    # pendulum mass
    l = 0.8    # pendulum length
    g = 9.81

    # dynamics: standard nonlinear pendulum-on-cart
    denom = M + m - m*cos(theta)**2

    x_dot     = v
    v_dot     = (F + m*sin(theta)*(l*omega**2 + g*cos(theta))) / denom
    theta_dot = omega
    omega_dot = (-F*cos(theta) - m*l*omega**2*cos(theta)*sin(theta)
                 - (M+m)*g*sin(theta)) / (l*denom)

    f_expl = vertcat(x_dot, v_dot, theta_dot, omega_dot)

    model = AcadosModel()
    model.name = model_name
    model.x = x_sym
    model.xdot = SX.sym('xdot', 4)
    model.u = u_sym
    model.z = SX.sym('z', 0, 0)
    model.p = SX.sym('p', 0, 0)
    model.f_expl_expr = f_expl
    # f_impl 안 쓰니까 설정 안 해도 됨

    return model


def diag_vec(diag_vals):
    """CasADi DM을 안 쓰고 numpy로만 대각 행렬 생성."""
    n = len(diag_vals)
    M = np.zeros((n, n))
    for i, v in enumerate(diag_vals):
        M[i, i] = v
    return M


def create_ocp():
    ocp = AcadosOcp()
    model = export_pendulum_model()
    ocp.model = model

    # dimensions
    nx = 4
    nu = 1
    N  = 40          # horizon steps
    Tf = 2.0         # horizon length [s]

    # ✨ 새 API: N은 options.N_horizon 에 넣기
    ocp.solver_options.N_horizon = N
    ocp.solver_options.tf = Tf

    # cost: LINEAR_LS
    ocp.cost.cost_type   = "LINEAR_LS"
    ocp.cost.cost_type_e = "LINEAR_LS"

    Q = diag_vec([10.0, 1.0, 50.0, 1.0])
    R = diag_vec([0.1])

    ny  = nx + nu
    ny_e = nx

    # Vx, Vu, Vx_e
    Vx   = np.zeros((ny, nx))
    Vu   = np.zeros((ny, nu))
    Vx_e = np.zeros((ny_e, nx))

    Vx[:nx, :nx] = np.eye(nx)
    Vu[nx:, 0]   = 1.0
    Vx_e[:, :]   = np.eye(nx)

    ocp.cost.Vx   = Vx
    ocp.cost.Vu   = Vu
    ocp.cost.Vx_e = Vx_e

    # W, W_e
    W   = np.zeros((ny, ny))
    W_e = np.zeros((ny_e, ny_e))

    for i in range(nx):
        W[i, i]   = Q[i, i]
        W_e[i, i] = Q[i, i]
    W[nx, nx] = R[0, 0]

    ocp.cost.W   = W
    ocp.cost.W_e = W_e

    ocp.cost.yref   = np.zeros((ny,))
    ocp.cost.yref_e = np.zeros((ny_e,))

    # input constraints
    Fmax = 15.0
    ocp.constraints.lbu = np.array([-Fmax])
    ocp.constraints.ubu = np.array([ Fmax])
    ocp.constraints.idxbu = np.array([0])

    # initial state
    ocp.constraints.x0 = np.array([0.0, 0.0, 0.2, 0.0])

    # solver options
    ocp.solver_options.qp_solver        = "PARTIAL_CONDENSING_HPIPM"
    ocp.solver_options.hessian_approx   = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type  = "ERK"
    ocp.solver_options.nlp_solver_type  = "SQP_RTI"
    ocp.solver_options.qp_solver_cond_N = N
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.sim_method_num_steps  = 1

    return ocp


def main():
    ocp = create_ocp()

    print("Generating solver + exporting C code...")

    # 현재 acados_template 버전에 맞는 방식: generate=True
    solver = AcadosOcpSolver(
        ocp,
        json_file="acados_ocp_pendulum.json",
        generate=True,
    )

    out_dir = os.path.join(os.getcwd(), "c_generated_code")
    print(f"✔ DONE. C code exported to: {out_dir}")


if __name__ == "__main__":
    main()
