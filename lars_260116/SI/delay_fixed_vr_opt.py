# =========================================================
# SINDy-like Sway–Yaw ID with Separate Fixed Delays (multi-trajectory)
# - Constrained LSQ:
#     vdot:  Yv<=0, Yvv<=0, Yr<=0,  a2>=0
#     rdot:  Nr<=0, Nrr<=0, Nv<=0,  a3>=0
# - Input term: (Tt^2) * sin(2.0 * Ts) in both v̇ and ṙ
# - Plots (one window per trajectory): v/r measured vs predicted + inputs
# =========================================================
import numpy as np
import matplotlib.pyplot as plt
import glob, random
from scipy.optimize import lsq_linear

# -----------------------------
# 0) Utils / Seed
# -----------------------------
def set_seed(seed=42):
    random.seed(seed); np.random.seed(seed)
set_seed(42)

# -----------------------------
# 0.1) Conversion functions
# -----------------------------

# -----------------------------
# 1) Load CSVs (v, r, throttle, steering)
# -----------------------------
def load_csv_data_vr(glob_pat):
    files = sorted(glob.glob(glob_pat))
    if len(files) == 0:
        raise RuntimeError(f"No CSV files found at {glob_pat}")
    Xvr_list, Tt_list, Ts_list = [], [], []
    for f in files:
        data = np.genfromtxt(f, delimiter=',', names=True)
        need = {'v','r','applied_thrust','applied_steer_rad'}
        if not need.issubset(set(data.dtype.names)):
            raise RuntimeError(f"{f} must contain columns: {need}")
        v = np.asarray(data['v'], dtype=float).reshape(-1,1)
        r = np.asarray(data['r'], dtype=float).reshape(-1,1)
        Xvr = np.hstack([v, r])
        tau_T = np.array([x/100.0 for x in data['applied_thrust']], dtype=float)
        tau_S = np.array([x for x in data['applied_steer_rad']], dtype=float)
        Xvr_list.append(Xvr)
        Tt_list.append(tau_T)
        Ts_list.append(tau_S)
    return Xvr_list, Tt_list, Ts_list

# -----------------------------
# 2) Align with separate delays (throttle vs steering)
#    States start at t = max(delay_T, delay_S)
#    Tt(t) -> Tt(t - delay_T), Ts(t) -> Ts(t - delay_S)
# -----------------------------
def align_separate_delay_vr(Xvr_list, Tt_list, Ts_list, delay_T, delay_S):
    Xs, Us = [], []
    dropped = 0
    max_delay = max(delay_T, delay_S)
    for Xvr, Tt, Ts in zip(Xvr_list, Tt_list, Ts_list):
        N = len(Xvr)
        if N <= max_delay:
            dropped += 1
            continue
        X_aligned = Xvr[max_delay:]  # drop first max_delay states

        # Align Tt and Ts windows separately to match the state times
        Tt_aligned = Tt[max_delay - delay_T : N - delay_T]
        Ts_aligned = Ts[max_delay - delay_S : N - delay_S]

        M = min(len(X_aligned), len(Tt_aligned), len(Ts_aligned))
        if M <= 1:
            dropped += 1
            continue

        Xs.append(X_aligned[:M])  # shape (M, 2) -> [v, r]
        Us.append(np.vstack([Tt_aligned[:M], Ts_aligned[:M]]).T)  # shape (M, 2) -> [Tt, Ts]
    if len(Xs) == 0:
        raise ValueError("All trajectories became empty after aligning for delays.")
    if dropped > 0:
        print(f"[WARN] Dropped {dropped} trajectory(ies) shorter than delays.")
    return Xs, Us, max_delay

# -----------------------------
# 3) Smoothed finite difference (per column)
# -----------------------------
def smoothed_finite_difference_col(x, dt, window=5):
    x = np.asarray(x).reshape(-1)
    if window > 1:
        k = window
        pad = np.pad(x, (k//2, k-1-k//2), mode='edge')
        xs = np.convolve(pad, np.ones(k)/k, mode='valid')
    else:
        xs = x
    xdot = np.zeros_like(xs)
    xdot[1:-1] = (xs[2:] - xs[:-2])/(2*dt)
    xdot[0]    = (xs[1] - xs[0])/dt
    xdot[-1]   = (xs[-1] - xs[-2])/dt
    return xdot

# -----------------------------
# 4) Build Θ and ẋ for v and r
#    vdot features: [v, |v|v, r, (Tt^2) sin(2Ts)]
#    rdot features: [r, |r|r, v, (Tt^2) sin(2Ts)]
# -----------------------------
def build_regression_matrices_vr(Xs, Us, dt, smooth_window=5, sin_factor=1.0):
    Theta_v_list, vdot_list = [], []
    Theta_r_list, rdot_list = [], []
    for X, U in zip(Xs, Us):
        v = X[:,0]; r = X[:,1]
        Tt, Ts = U[:,0], U[:,1]
        u_term = (Tt**2) * np.sin(Ts)

        # vdot model
        f_v = np.column_stack([
            v,                 # Yv
            np.abs(v)*v,       # Yvv
            r,                 # Yr
            u_term             # a2
        ])
        vdot = smoothed_finite_difference_col(v, dt, window=smooth_window)

        # rdot model
        f_r = np.column_stack([
            r,                 # Nr
            np.abs(r)*r,       # Nrr
            v,                 # Nv
            u_term             # a3
        ])
        rdot = smoothed_finite_difference_col(r, dt, window=smooth_window)

        m = min(len(vdot), len(f_v), len(rdot), len(f_r))
        Theta_v_list.append(f_v[:m,:])
        vdot_list.append(vdot[:m])
        Theta_r_list.append(f_r[:m,:])
        rdot_list.append(rdot[:m])

    A_v = np.vstack(Theta_v_list)
    b_v = np.concatenate(vdot_list)
    A_r = np.vstack(Theta_r_list)
    b_r = np.concatenate(rdot_list)
    return A_v, b_v, A_r, b_r

# -----------------------------
# 5) Constrained LSQ for (v̇, ṙ), with tiny ridge regularization
# -----------------------------
def fit_constrained_vr(Xs, Us, dt, l2=1e-8, smooth_window=5, sin_factor=2.0):
    A_v, b_v, A_r, b_r = build_regression_matrices_vr(
        Xs, Us, dt, smooth_window=smooth_window, sin_factor=sin_factor
    )

    def solve_block(A, b, lb, ub, l2):
        # scale columns
        col_scale = np.maximum(A.std(axis=0), 1e-8)
        A_s = A / col_scale

        # optional ridge via augmentation
        if l2 > 0:
            lam = np.sqrt(l2)
            A_aug = np.vstack([A_s, lam*np.eye(A_s.shape[1])])
            b_aug = np.concatenate([b, np.zeros(A_s.shape[1])])
        else:
            A_aug, b_aug = A_s, b

        res = lsq_linear(A_aug, b_aug, bounds=(lb, ub), method='trf', verbose=0)
        coef = res.x / col_scale
        cost = res.cost / len(b)
        return coef, cost, res

    # vdot bounds: [Yv, Yvv, Yr, a2] -> (<=0, <=0, <=0, >=0)
    lb_v = np.array([-np.inf, -np.inf, -np.inf, 0.0])
    ub_v = np.array([ 0.0,     0.0,     0.0,    np.inf])
    coef_v, cost_v, _ = solve_block(A_v, b_v, lb_v, ub_v, l2)

    # rdot bounds: [Nr, Nrr, Nv, a3] -> (<=0, <=0, <=0, >=0)
    lb_r = np.array([-np.inf, -np.inf, -np.inf, -np.inf])
    ub_r = np.array([ 0.0,     0.0,     0.0,    0.0])
    coef_r, cost_r, _ = solve_block(A_r, b_r, lb_r, ub_r, l2)

    (Yv, Yvv, Yr, a2) = coef_v.tolist()
    (Nr, Nrr, Nv, a3) = coef_r.tolist()
    return (Yv, Yvv, Yr, a2, Nr, Nrr, Nv, a3), (cost_v, cost_r)

# -----------------------------
# 6) Residual score (MSE on v̇ and ṙ)
# -----------------------------
def residual_score_vr(params, Xs, Us, dt, smooth_window=5, sin_factor=2.0):
    Yv, Yvv, Yr, a2, Nr, Nrr, Nv, a3 = params
    errs_v, errs_r = [], []
    for X, U in zip(Xs, Us):
        v = X[:,0]; r = X[:,1]
        Tt, Ts = U[:,0], U[:,1]
        u_term = (Tt**2) * np.sin(sin_factor * Ts)

        vdot = smoothed_finite_difference_col(v, dt, window=smooth_window)
        rdot = smoothed_finite_difference_col(r, dt, window=smooth_window)

        vdot_hat = Yv*v + Yvv*np.abs(v)*v + Yr*r + a2*u_term
        rdot_hat = Nr*r + Nrr*np.abs(r)*r + Nv*v + a3*u_term

        m = min(len(vdot), len(vdot_hat), len(rdot), len(rdot_hat))
        errs_v.append(np.mean((vdot[:m] - vdot_hat[:m])**2))
        errs_r.append(np.mean((rdot[:m] - rdot_hat[:m])**2))
    return float(np.mean(errs_v)) if errs_v else np.inf, float(np.mean(errs_r)) if errs_r else np.inf

# -----------------------------
# 7) Coupled rollout simulation (Euler, no delays inside)
# -----------------------------
def simulate_vr(v0, r0, U_eff, dt, params, sin_factor=2.0):
    Yv, Yvv, Yr, a2, Nr, Nrr, Nv, a3 = params
    N = len(U_eff)
    t = np.arange(N)*dt
    v = np.zeros(N); r = np.zeros(N)
    v[0] = float(v0); r[0] = float(r0)
    Tt, Ts = U_eff[:,0], U_eff[:,1]
    for k in range(1, N):
        u_term = (Tt[k-1]**2) * np.sin(sin_factor * Ts[k-1])
        vdot = Yv*v[k-1] + Yvv*np.abs(v[k-1])*v[k-1] + Yr*r[k-1] + a2*u_term
        rdot = Nr*r[k-1] + Nrr*np.abs(r[k-1])*r[k-1] + Nv*v[k-1] + a3*u_term
        v[k] = v[k-1] + dt*vdot
        r[k] = r[k-1] + dt*rdot
    return t, v.reshape(-1,1), r.reshape(-1,1)

# -----------------------------
# 8) Main
# -----------------------------
if __name__ == "__main__":
    # CSV_GLOB = "dataset/vr/*.csv"   # <-- change path as needed
    # CSV_GLOB = "check/*.csv"
    CSV_GLOB =  "/home/user/aura_ws/src/lars/SI/sway_yaw/*.csv"
    dt = 0.1

    # Separate fixed delays (seconds -> steps)
    delay_T_sec = 0.1
    delay_S_sec = 0.1
    delay_T = int(round(delay_T_sec / dt))   # throttle delay
    delay_S = int(round(delay_S_sec / dt))   # steering delay

    # Load & align
    Xvr_list, Tt_list, Ts_list = load_csv_data_vr(CSV_GLOB)
    Xs, Us, max_delay = align_separate_delay_vr(Xvr_list, Tt_list, Ts_list, delay_T, delay_S)

    # Fit constrained LSQ
    params, (cost_v, cost_r) = fit_constrained_vr(Xs, Us, dt, l2=1e-8, smooth_window=5, sin_factor=2.0)
    Yv, Yvv, Yr, a2, Nr, Nrr, Nv, a3 = params
    v_mse, r_mse = residual_score_vr(params, Xs, Us, dt, smooth_window=5, sin_factor=2.0)

    print("\n=== SWAY–YAW CONSTRAINED FIT (T delay=%.1fs, S delay=%.1fs) ===" % (delay_T_sec, delay_S_sec))
    print(f"vdot:  Yv={Yv:.6f} (<=0), Yvv={Yvv:.6f} (<=0), Yr={Yr:.6f} (<=0), a2={a2:.6f} (>=0)")
    print(f"rdot:  Nr={Nr:.6f} (<=0), Nrr={Nrr:.6f} (<=0), Nv={Nv:.6f} (<=0), a3={a3:.6f} (>=0)")
    print(f"LSQ costs (per-sample): v: {cost_v:.6e} | r: {cost_r:.6e}")
    print(f"Residual MSE: v̇: {v_mse:.6e} | ṙ: {r_mse:.6e}")
    print("Models:")
    print("  vdot = Yv*v + Yvv*|v|v + Yr*r + a2*(Tt^2*sin(2.0*Ts))")
    print("  rdot = Nr*r + Nrr*|r|r + Nv*v + a3*(Tt^2*sin(2.0*Ts))")

    # -------------------------
    # Plots (per trajectory): states and inputs in the SAME window
    # -------------------------
    for idx, (X_raw, Tt_raw, Ts_raw) in enumerate(zip(Xvr_list, Tt_list, Ts_list)):
        N = len(X_raw)
        if N <= max_delay:
            continue

        # Effective aligned segment for this trajectory:
        X_eff  = X_raw[max_delay:]  # (Nv, 2) -> [v, r]
        Tt_eff = Tt_raw[max_delay - delay_T : N - delay_T]
        Ts_eff = Ts_raw[max_delay - delay_S : N - delay_S]
        M = min(len(X_eff), len(Tt_eff), len(Ts_eff))
        if M <= 1: 
            continue
        X_eff  = X_eff[:M]
        U_eff  = np.vstack([Tt_eff[:M], Ts_eff[:M]]).T
        t_eff  = np.arange(M)*dt + max_delay*dt

        # Simulate coupled (v,r) with identified parameters
        t_sim, v_sim, r_sim = simulate_vr(X_eff[0,0], X_eff[0,1], U_eff, dt, params, sin_factor=2.0)

        # Build delayed control traces for visualization
        t_full = np.arange(N)*dt
        Tt_delayed_plot = np.concatenate([np.full(delay_T, np.nan), Tt_raw[:-delay_T]])
        Ts_delayed_plot = np.concatenate([np.full(delay_S, np.nan), Ts_raw[:-delay_S]])

        # --- One figure with 3 rows: v, r, inputs ---
        fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=False)
        fig.suptitle(f"Trajectory {idx+1}: Throttle delay={delay_T_sec:.1f}s, Steering delay={delay_S_sec:.1f}s")

        # Row 1: sway speed v
        axes[0].plot(t_eff, X_eff[:,0], 'k', lw=2, label='Measured v (aligned)')
        axes[0].plot(t_eff, v_sim[:,0], 'r--', lw=2, label='Predicted v')
        axes[0].axvline(max_delay*dt, color='gray', ls='--', label='start of aligned window')
        axes[0].set_ylabel('Sway v')
        axes[0].grid(True); axes[0].legend()

        # Row 2: yaw rate r
        axes[1].plot(t_eff, X_eff[:,1], 'k', lw=2, label='Measured r (aligned)')
        axes[1].plot(t_eff, r_sim[:,0], 'r--', lw=2, label='Predicted r')
        axes[1].axvline(max_delay*dt, color='gray', ls='--')
        axes[1].set_ylabel('Yaw r')
        axes[1].grid(True); axes[1].legend()

        # Row 3: inputs (original vs delayed)
        axes[2].plot(t_full, Tt_raw, 'b', alpha=0.35, label='Throttle (original)')
        axes[2].plot(t_full, Tt_delayed_plot, 'b', lw=2, label='Throttle (delayed)')
        axes[2].plot(t_full, Ts_raw, 'g', alpha=0.35, label='Steering (original)')
        axes[2].plot(t_full, Ts_delayed_plot, 'g', lw=2, label='Steering (delayed)')
        axes[2].axvline(max_delay*dt, color='gray', ls='--')
        axes[2].set_xlabel('Time [s]'); axes[2].set_ylabel('Inputs (scaled)')
        axes[2].grid(True); axes[2].legend(loc='best')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.show()
