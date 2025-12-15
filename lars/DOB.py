import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64MultiArray
import math
import time
import numpy as np


def DOB(F,state, state_estim, param_filtered, param_estim, dt):


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

    w_cutoff = 1.0
    gain = -1.0


    # set up states & controls
    u    = state[3]
    v    = state[4]
    r    = state[5]
    delta  = state[6]
    # F_eff  = state[7]
    F_eff  = F
    eps = 1e-6

    # print(u, v, r, F)
 
    # # F = F_cmd 
    # Deadzone
    s = 25
    k = 8
    a1 = 2.2*2.2
    a2 = 2.2*2.2
    b11 = 1.0
    b22 = 1.0
    T = ((1/(1+np.exp(s*F_eff)))*(b11*F_eff + np.tanh(k*F_eff)*a1) + (1/(1+np.exp(-s*F_eff)))*(b22*F_eff + np.tanh(k*F_eff)*a2))

    
    f_usv = np.array([( - Xu*u - Xuu * np.sqrt(u * u + eps) * u + alpha1*T*np.cos(alpha2*delta)) ,
                    ( -Yv*v - Yr*r - Yvv * np.sqrt(v * v + eps) * v + alpha3*T*np.sin(alpha2*delta)) ,
                    ( - Nr*r - Nv*v - Nrr * np.sqrt(r * r + eps) * r - alpha4*T*np.sin(alpha2*delta))
                    ])
    # print(f_usv[0])

    uvr = np.array([u, v, r])

    x_error = state_estim - uvr 
    
    
    xdot = param_estim + f_usv + gain*x_error

    x_t_plus = xdot*dt + state_estim
    


    pi = (1/gain)*(np.exp(gain*dt)-1.0)
    param_estim = -np.exp(gain*dt)*x_error/pi
    
    before_param_filtered = param_filtered
    param_filtered = before_param_filtered*math.exp(-w_cutoff*dt) - param_estim*(1-math.exp(-w_cutoff*dt))

    
    return x_t_plus, param_estim, param_filtered


