# =========================================================
# SINDy-like Surge ID with Separate Fixed Delays (multi-trajectory)
# - True constrained LSQ: c0<=0 (Xu), c1<=0 (Xuu), c2>=0 (alpha1)
# - Throttle delay = 3.0 s, Steering delay = 2.0 s
# - Plots measured vs predicted surge & original vs delayed inputs
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
def convertThrustToPwm(thrust):
    if thrust <= 1500:
        pwm1 = (thrust - 1450) * 0.26
    else:
        pwm1 = (thrust - 1550) * 0.26
    return max(pwm1, 1.3)

def convertSteeringToPwm(pwm):
    if pwm >= 2000:
        steer = 300
    elif 1550 < pwm < 2000:
        steer = (pwm - 1550) / 1.6667
    elif 1000 <= pwm <= 1450:
        steer = (pwm - 1450) / 1.6667
    elif pwm <= 1000:
        steer = -300
    else:
        steer = 0
    return steer

# -----------------------------
# 1) Load CSVs (u, throttle, steering)
# -----------------------------
def load_csv_data(glob_pat):
    files = sorted(glob.glob(glob_pat))
    if len(files) == 0:
        raise RuntimeError(f"No CSV files found at {glob_pat}")
    X_list, Tt_list, Ts_list = [], [], []
    for f in files:
        data = np.genfromtxt(f, delimiter=',', names=True)
        need = {'u','applied_thrust','applied_steer_rad'}
        if not need.issubset(set(data.dtype.names)):
            raise RuntimeError(f"{f} must contain columns: {need}")
        u = np.asarray(data['u'], dtype=float).reshape(-1,1)
        tau_T = np.array([x/100.0 for x in data['applied_thrust']], dtype=float)
        tau_S = np.array([x for x in data['applied_steer_rad']], dtype=float)
        X_list.append(u)
        Tt_list.append(tau_T)
        Ts_list.append(tau_S)
    return X_list, Tt_list, Ts_list

# -----------------------------
# 2) Align with separate delays
# -----------------------------
def align_separate_delay(X_list, Tt_list, Ts_list, delay_T, delay_S):
    Xs, Us = [], []
    dropped = 0
    max_delay = max(delay_T, delay_S)
    for u, Tt, Ts in zip(X_list, Tt_list, Ts_list):
        N = len(u)
        if N <= max_delay:
            dropped += 1
            continue
        X_aligned = u[max_delay:]
        Tt_aligned = Tt[max_delay - delay_T : N - delay_T]
        Ts_aligned = Ts[max_delay - delay_S : N - delay_S]
        M = min(len(X_aligned), len(Tt_aligned), len(Ts_aligned))
        if M <= 1:
            dropped += 1
            continue
        Xs.append(X_aligned[:M])
        Us.append(np.vstack([Tt_aligned[:M], Ts_aligned[:M]]).T)
    if len(Xs) == 0:
        raise ValueError("All trajectories became empty after aligning for delays.")
    if dropped > 0:
        print(f"[WARN] Dropped {dropped} trajectory(ies) shorter than delays.")
    return Xs, Us, max_delay

# -----------------------------
# 3) Finite difference (smoothed)
# -----------------------------
def smoothed_finite_difference(x, dt, window=5):
    # simple Savitzky-Golay-like smoothing via moving average, then finite diff
    x = np.asarray(x).reshape(-1,1)
    if window > 1:
        k = window
        pad = np.pad(x, ((k//2, k-1-k//2),(0,0)), mode='edge')
        x_s = np.convolve(pad[:,0], np.ones(k)/k, mode='valid').reshape(-1,1)
    else:
        x_s = x
    xdot = np.zeros_like(x_s)
    xdot[1:-1,0] = (x_s[2:,0] - x_s[:-2,0])/(2*dt)
    xdot[0,0]    = (x_s[1,0] - x_s[0,0])/dt
    xdot[-1,0]   = (x_s[-1,0] - x_s[-2,0])/dt
    return xdot

# -----------------------------
# 4) Build features Θ and targets ẋ for all trajectories
#    Θ columns: [u, |u|u, Tt^2 cos(2.0 Ts)]
# -----------------------------
def build_regression_matrices(Xs, Us, dt, smooth_window=5):
    Theta_list, xdot_list = [], []
    for X, U in zip(Xs, Us):
        u = X[:,0]
        Tt = U[:,0]; Ts = U[:,1]
        f1 = u
        f2 = np.abs(u)*u
        f3 = (Tt**2)*np.cos(Ts)
        Theta = np.column_stack([f1, f2, f3])
        xdot = smoothed_finite_difference(u, dt, window=smooth_window)[:,0]
        # align lengths (safe)
        m = min(len(xdot), len(Theta))
        Theta_list.append(Theta[:m,:])
        xdot_list.append(xdot[:m])
    A = np.vstack(Theta_list)
    b = np.concatenate(xdot_list)
    return A, b

# -----------------------------
# 5) Constrained LSQ fit (with optional tiny ridge)
# -----------------------------
def fit_constrained(Xs, Us, dt, l2=1e-8, smooth_window=5):
    # Build stacked Θ, ẋ
    A, b = build_regression_matrices(Xs, Us, dt, smooth_window=smooth_window)

    # Column scaling for numerical conditioning
    col_scale = np.maximum(A.std(axis=0), 1e-8)
    A_s = A / col_scale

    # Optional ridge via row augmentation: sqrt(l2)*I
    if l2 > 0:
        lam = np.sqrt(l2)
        A_aug = np.vstack([A_s, lam*np.eye(A_s.shape[1])])
        b_aug = np.concatenate([b, np.zeros(A_s.shape[1])])
    else:
        A_aug, b_aug = A_s, b

    # Bounds: c0<=0, c1<=0, c2>=0  -> lb=[-inf,-inf,0], ub=[0,0,inf]
    lb = np.array([-np.inf, -np.inf, 0.0])
    ub = np.array([ 0.0,     0.0,    np.inf])

    res = lsq_linear(A_aug, b_aug, bounds=(lb, ub), method='trf', verbose=0)
    coef_scaled = res.x
    coef = coef_scaled / col_scale  # undo scaling
    Xu, Xuu, alpha1 = coef.tolist()
    return (Xu, Xuu, alpha1), res.cost/len(b)

# -----------------------------
# 6) Residual score (per-trajectory MSE on ẋ)
# -----------------------------
def residual_score_params(params, Xs, Us, dt, smooth_window=5):
    Xu, Xuu, alpha1 = params
    errs = []
    for X, U in zip(Xs, Us):
        u = X[:,0]
        Tt, Ts = U[:,0], U[:,1]
        xdot = smoothed_finite_difference(u, dt, window=smooth_window)[:,0]
        f1 = u
        f2 = np.abs(u)*u
        f3 = (Tt**2)*np.cos(2.0*Ts)
        xdot_hat = Xu*f1 + Xuu*f2 + alpha1*f3
        m = min(len(xdot), len(xdot_hat))
        errs.append(np.mean((xdot[:m] - xdot_hat[:m])**2))
    return float(np.mean(errs)) if errs else np.inf

# -----------------------------
# 7) Rollout simulation
# -----------------------------
def simulate_with_identified(u0, U_eff, dt, params):
    Xu, Xuu, alpha1 = params
    N = len(U_eff)
    t = np.arange(N)*dt
    u = np.zeros((N,1)); u[0,0] = float(u0)
    Tt, Ts = U_eff[:,0], U_eff[:,1]
    for k in range(1, N):
        uk = u[k-1,0]
        du = Xu*uk + Xuu*np.abs(uk)*uk + alpha1*(Tt[k-1]**2)*np.cos(2.0*Ts[k-1])
        u[k,0] = uk + dt*du
    return t, u

# -----------------------------
# 8) Main
# -----------------------------
if __name__ == "__main__":
    CSV_GLOB = "/home/user/aura_ws/src/lars/SI/surge/*.csv"
    dt = 0.1
    delay_T = int(round(0.1 / dt))   # throttle delay = 3.0 s
    delay_S = int(round(1.0 / dt))   # steering delay = 2.0 s

    # Load & align
    X_list, Tt_list, Ts_list = load_csv_data(CSV_GLOB)
    Xs, Us, max_delay = align_separate_delay(X_list, Tt_list, Ts_list, delay_T, delay_S)

    # Fit constrained LSQ
    params, global_cost = fit_constrained(Xs, Us, dt, l2=1e-8, smooth_window=5)
    Xu_hat, Xuu_hat, alpha1_hat = params
    rescore = residual_score_params(params, Xs, Us, dt, smooth_window=5)

    print("\n=== FIXED DELAY CONSTRAINED FIT (Throttle=3.0s, Steering=2.0s) ===")
    print(f"Xu      : {Xu_hat:.6f}  (<= 0)")
    print(f"Xuu     : {Xuu_hat:.6f} (<= 0)")
    print(f"alpha1  : {alpha1_hat:.6f} (>= 0)")
    print(f"Global LSQ cost (per-sample): {global_cost:.6e}")
    print(f"Residual ẋ MSE: {rescore:.6e}\n")
    print("Identified model: du = Xu*u + Xuu*|u|u + alpha1*(Tt^2*cos(2.0*Ts))")

    # -------------------------
    # Plots (per trajectory)
    # -------------------------
    for idx, (u_raw, Tt_raw, Ts_raw) in enumerate(zip(X_list, Tt_list, Ts_list)):
        N = len(u_raw)
        if N <= max_delay:
            continue

        # Effective aligned segment
        X_eff  = u_raw[max_delay:]
        Tt_eff = Tt_raw[max_delay - delay_T : N - delay_T]
        Ts_eff = Ts_raw[max_delay - delay_S : N - delay_S]
        M = min(len(X_eff), len(Tt_eff), len(Ts_eff))
        if M <= 1: continue
        X_eff = X_eff[:M]
        U_eff = np.vstack([Tt_eff[:M], Ts_eff[:M]]).T

        # Predict rollout
        t_eff, u_sim = simulate_with_identified(X_eff[0,0], U_eff, dt, params)

        # Build delayed-control traces for visualization
        t_full = np.arange(N)*dt
        Tt_delayed_plot = np.concatenate([np.full(delay_T, np.nan), Tt_raw[:-delay_T]])
        Ts_delayed_plot = np.concatenate([np.full(delay_S, np.nan), Ts_raw[:-delay_S]])

        # --- Plot surge ---
        plt.figure(figsize=(12,6))
        plt.subplot(2,1,1)
        plt.plot(t_eff + max_delay*dt, X_eff[:,0], "k", lw=2, label="Measured u (aligned)")
        plt.plot(t_eff + max_delay*dt, u_sim[:,0], "r--", lw=2, label="Predicted u")
        plt.axvline(max_delay*dt, color="gray", ls="--", label="start of aligned window")
        plt.ylabel("Surge speed u")
        plt.title(f"Trajectory {idx+1}: Throttle=3.0s, Steering=2.0s")
        plt.legend(); plt.grid(True)

        # --- Plot inputs ---
        plt.subplot(2,1,2)
        plt.plot(t_full, Tt_raw, "b", alpha=0.35, label="Throttle (original)")
        plt.plot(t_full, Tt_delayed_plot, "b", lw=2, label="Throttle (delayed)")
        plt.plot(t_full, Ts_raw, "g", alpha=0.35, label="Steering (original)")
        plt.plot(t_full, Ts_delayed_plot, "g", lw=2, label="Steering (delayed)")
        plt.axvline(max_delay*dt, color="gray", ls="--")
        plt.xlabel("Time [s]"); plt.ylabel("Inputs (scaled)")
        plt.legend(loc="best"); plt.grid(True)
        plt.tight_layout()

    plt.show()
