import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Load log
# =========================
CSV_PATH = "/home/user/aura_ws/src/lars/SI/candidate/13.csv"
df = pd.read_csv(CSV_PATH)

# time (relative, seconds)
t = df["t_sec"].to_numpy()
t = t - t[0]

# =========================
# 1) (x, y) trajectory
# =========================
x = df["x"].to_numpy()
y = df["y"].to_numpy()

plt.figure()
plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("y")
plt.title("Trajectory (x-y)")
plt.axis("equal")
plt.grid(True)

# =========================
# 2) u, v, r history
# =========================
plt.figure()
plt.plot(t, df["u"], label="u")
plt.xlabel("t [s]")
plt.ylabel("u")
plt.title("u history")
plt.grid(True)
plt.legend()

plt.figure()
plt.plot(t, df["v"], label="v")
plt.xlabel("t [s]")
plt.ylabel("v")
plt.title("v history")
plt.grid(True)
plt.legend()

plt.figure()
plt.plot(t, df["r"], label="r")
plt.xlabel("t [s]")
plt.ylabel("r")
plt.title("r history")
plt.grid(True)
plt.legend()

# =========================
# 3) Control input history
#   - plots both applied and target (if present)
# =========================
# Steering (UI units + rad if available)
plt.figure()
if "applied_steer_ui" in df.columns:
    plt.plot(t, df["applied_steer_ui"], label="applied_steer_ui")
if "target_steer_ui" in df.columns:
    plt.plot(t, df["target_steer_ui"], "--", label="target_steer_ui")
plt.xlabel("t [s]")
plt.ylabel("steer (UI)")
plt.title("Steering input history (UI)")
plt.grid(True)
plt.legend()

if "applied_steer_rad" in df.columns:
    plt.figure()
    plt.plot(t, df["applied_steer_rad"], label="applied_steer_rad")
    plt.xlabel("t [s]")
    plt.ylabel("steer [rad]")
    plt.title("Steering input history (rad)")
    plt.grid(True)
    plt.legend()

# Thrust (applied + target)
plt.figure()
if "applied_thrust" in df.columns:
    plt.plot(t, df["applied_thrust"], label="applied_thrust")
if "target_thrust" in df.columns:
    plt.plot(t, df["target_thrust"], "--", label="target_thrust")
plt.xlabel("t [s]")
plt.ylabel("thrust")
plt.title("Thrust input history")
plt.grid(True)
plt.legend()

plt.show()
