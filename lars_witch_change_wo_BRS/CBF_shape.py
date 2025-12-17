import numpy as np
import matplotlib.pyplot as plt

def world_to_local(x, y, CD_x, CD_y, CD_psi):
    dx = x - CD_x
    dy = y - CD_y
    x_l =  np.cos(CD_psi)*dx + np.sin(CD_psi)*dy
    y_l = -np.sin(CD_psi)*dx + np.cos(CD_psi)*dy
    return x_l, y_l


def hopper_upper_cbf(x, y, CD_x, CD_y, CD_psi,
                     a, b, x0, y0):

    x_l, y_l = world_to_local(x, y, CD_x, CD_y, CD_psi)

    # 방향 반전 핵심
    y_upper = a * np.tanh(b * (x0 - x_l)) + y0
    return y_upper - y_l


def hopper_lower_cbf(x, y, CD_x, CD_y, CD_psi,
                     a, b, x0, y0):

    x_l, y_l = world_to_local(x, y, CD_x, CD_y, CD_psi)

    # 방향 반전 핵심
    y_lower = -a * np.tanh(b * (x0 - x_l)) - y0
    return y_l - y_lower



X, Y = np.meshgrid(
    np.linspace(-10, 10, 600),
    np.linspace(-8, 8, 600)
)

CD_x, CD_y = 0.0, 0.0
CD_psi = np.deg2rad(0.0)

a = 3.0     # hopper half-height
b = 0.3     # 곡률 (얼마나 급하게 열리는지)
x0 = -10.0    # hopper 목 위치
y0 = 4.0    # 중심 오프셋


hU = hopper_upper_cbf(X, Y, CD_x, CD_y, CD_psi, a, b, x0, y0)
hL = hopper_lower_cbf(X, Y, CD_x, CD_y, CD_psi, a, b, x0, y0)

feasible = np.minimum(hU, hL)

plt.figure(figsize=(6,6))
levels = np.linspace(-4, 4, 80)
cf = plt.contourf(X, Y, feasible, levels=levels, cmap='coolwarm')
plt.colorbar(cf, label='min(h_upper, h_lower)')

plt.contour(X, Y, hU, levels=[0], colors='k', linewidths=2)
plt.contour(X, Y, hL, levels=[0], colors='k', linewidths=2)

plt.plot(CD_x, CD_y, 'ko')
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.title('tanh 기반 Hopper CBF (두 개의 경계 곡선)')
plt.show()
