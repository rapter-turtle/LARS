import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64MultiArray, Float64
from aura_msg.msg import Waypoint, Parameter
from sensor_msgs.msg import Imu, NavSatFix
from math import atan2, sqrt, pi, cos, sin, log
from geographiclib.geodesic import Geodesic
from pyproj import Transformer
import casadi as ca
import numpy as np


def clamp(value, low, high):
    return max(low, min(value, high))

class ActuatorPublisher(Node):
    def __init__(self):
        super().__init__('actuator_publisher')

        # Publishers
        self.publisher_ = self.create_publisher(Float64MultiArray, '/actuator_outputs', 10)
        self.utm_publisher_ = self.create_publisher(Float64MultiArray, '/cbf', 10)

        # Subscribers
        self.create_subscription(Float64MultiArray, '/ekf/estimated_state', self.state_callback, 10)
        self.create_subscription(Float64, '/desired_velocity', self.velocity_callback, 10)

        # Timer
        self.timer_ = self.create_timer(0.1, self.timer_callback)


        x_actual_min = 289577.66
        x_actual_max = 291591.05
        y_actual_min = 4117065.30
        y_actual_max = 4118523.52
        
        station_keeping_point = np.array([(x_actual_min+x_actual_max)*0.5, (y_actual_min+y_actual_max)*0.5])        
        # Control parameters and state
        self.x = self.y = self.psi = self.u = self.v = self.r = 0.0
        self.x = station_keeping_point[0]
        self.y = station_keeping_point[1]
        self.desired_velocity = 0.0
        self.max_steer = 300
        self.max_thrust = 100
        self.max_thrust_diff = 100.0
        self.max_steer_diff = 100.0
        self.before_velocity_e = 0.0
        self.last_steering = 0.0
        self.last_thrust = 0.0
        self.gps_received = False
        self.imu_received = False

        # System parameter
        self.Xu = 0.10531
        self.Xuu = 0.018405
        self.Yv = 8.7185e-09
        self.Yvv = 0.39199
        self.Yn = 1.1508e-08
        self.Nr = 9.1124e-08
        self.Nrr = 5.2726
        self.Nv = 6.1558e-09
        self.b1 = 0.00058466
        self.b2 = 0.0040635
        self.b3 = 0.31094
        self.hp = 3.5
        self.bp = -3.5

        self.dt = 0.1
        self.eta_x_d = 0.0
        self.eta_y_d = 0.0
        self.a = 3.0
        self.b = 0.01

        self.alpha1 = 20.0
        self.alpha2 = 10.0
        self.alpha_clf = 1.0

        self.eta_ux = 0.0
        self.eta_uy = 0.0

        self.tau_x = 0.0
        self.tau_y = 0.0
        self.tau_u = 0.0
        self.tau_r = 0.0

        self.desired_u = 0.0
        self.desired_r = 0.0

        self.K_r = 10.0
        self.K_u = 20.0
        self.Kd_r = -10.0
        self.Kd_u = -20.0        
        self.eu_before = 0.0
        self.er_before = 0.0

        self.last_u = 0.0
        self.last_r = 0.0



        self.x0 = station_keeping_point[0] - 0.0
        self.y0 = station_keeping_point[1] + 50.0 

        self.dt = 0.1
        self.V0 = 2.5

    def state_callback(self, msg):
        self.x, self.y, self.psi, self.u, self.v, self.r = msg.data[:6]
        self.gps_received = True
        self.imu_received = True
        self.publish_utm_coordinates()

    def velocity_callback(self, msg):
        self.desired_velocity = msg.data / 1.94384  # knots to m/s

    def convert_steering_to_pwm(self, steer):
        if steer >= 300:
            return 2000.0
        elif steer >= 0:
            return 1500.0 + steer * 1.6667
        elif steer >= -300:
            return 1500.0 + steer * 1.6667
        else:
            return 1000.0

    def convert_thrust_to_pwm(self, thrust):
        if thrust <= 0:
            return 1500.0
        else:
            pwm = 3.9 * thrust + 1550.0
            return clamp(pwm, 1500.0, 2000.0)


    def timer_callback(self):

        self.y0 = self.y0 + self.V0*self.dt

        psi = self.psi
        psi = psi - 2*pi if psi > pi else psi
        psi = psi + 2*pi if psi < -pi else psi



        # CLF-CBF optimization
        opti = ca.Opti()

        uu = opti.variable()  # First control input (thrust)
        rr = opti.variable()  # Second control input (steering)
        s1 = opti.variable()  # Slack variable
        s2 = opti.variable()
        s3 = opti.variable()

        # s4 = opti.variable()
        # Objective function

        opti.set_initial(uu, 0.0)         # initialize thrust with current u
        opti.set_initial(rr, 0.0)         # initialize steering with current r
        opti.set_initial(s1, 0.0)
        opti.set_initial(s2, 0.0)            # initialize slack variables
        opti.set_initial(s3, 0.0)

        slack_weight = 10000000.0

        # opti.minimize((uu - self.u)**2 + 0.01*(rr - self.r)**2 + slack_weight*s2**2 + slack_weight*s3**2 )
        opti.minimize((uu - self.u)**2 + 100000*(rr - self.r)**2 + slack_weight*s1**2 + slack_weight*s2**2 + slack_weight*s3**2 )


        ## Front CBF
        etaVx = uu*ca.cos(psi) - self.v*ca.sin(psi) - rr*self.hp*ca.sin(psi)
        etaVy = uu*ca.sin(psi) + self.v*ca.cos(psi) + rr*self.hp*ca.cos(psi) 
        etax = self.x + self.hp*ca.cos(psi)
        etay = self.y + self.hp*ca.sin(psi)
        guide_cbf = -(etay - self.y0) - self.a*log((etax - self.x0)*(etax - self.x0) + 1)
        guide_cbf_dot = -(etaVy - self.V0) - 2*self.a*(etax - self.x0)*etaVx/((etax - self.x0)*(etax - self.x0) + 1)
  
        ## CLF
        clf_x = (etax - self.x0)*(etax - self.x0) 
        clf_y = (etay - self.y0)*(etay - self.y0)

        opti.subject_to(uu >= -0.1)
        opti.subject_to(uu <= 4.0)        
        opti.subject_to(rr >= -0.2)
        opti.subject_to(rr <= 0.2)    
        opti.subject_to(rr >= self.desired_r - 0.02)
        opti.subject_to(rr <= self.desired_r + 0.02)           
        opti.subject_to(self.alpha1*guide_cbf + guide_cbf_dot + 1000000+ s1 >= 0.0)     
        opti.subject_to(self.alpha_clf*clf_x + 2*(etax - self.x0)*etaVx + s2 <= 0.0)
        opti.subject_to(self.alpha_clf*clf_y + 2*(etay - self.y0)*(etaVy - self.V0) + s3 <= 0.0)

        # print(self.alpha_clf*clf_x + 2*(etax - self.x0)*(self.u*cos(psi) - self.v*sin(psi) - self.r*self.hp*sin(psi)))
        opti.solver('ipopt')
        sol = opti.solve()
        # optimal_u = opti.value(uu)

        self.desired_r = opti.value(rr)
        self.desired_u = opti.value(uu)
        eu = (self.desired_u - self.u)
        er = (-self.desired_r - self.r)

        self.last_u = self.desired_u
        self.last_r = self.desired_r

        Ku = 200.0
        Kr = 700.0
        
        proposed_thrust = Ku*eu + self.Kd_u*(eu - self.eu_before)
        steer_input = Kr*er + self.Kd_r*(er - self.er_before)


        steer_input = clamp(steer_input, -self.max_steer, self.max_steer)
        
        steer_change = clamp(steer_input - self.last_steering, -self.max_steer_diff, self.max_steer_diff)
        thrust_change = clamp(proposed_thrust - self.last_thrust, -self.max_thrust_diff, self.max_thrust_diff)

        steer = self.last_steering + steer_change
        thrust = self.last_thrust + thrust_change
        self.last_steering = steer
        self.last_thrust = thrust
        self.eu_before = (self.desired_u - self.u)
        self.er_before = (self.desired_r - self.r)

        # print("tu = ", self.b1*proposed_thrust*proposed_thrust*cos(self.b2*optimal_p))
        # print("tr = ", -self.b3*self.b1*proposed_thrust*proposed_thrust*sin(self.b2*optimal_p))
        print("u = ", self.u, "desired u = ", self.desired_u)
        print("r = ", self.r, "desired r = ", self.desired_r)
        # print("rudder angle = ", optimal_p)
        # print("s1 : ", allocation.value(alloc_s1))

        pwm_steer = self.convert_steering_to_pwm(steer)
        pwm_thrust = self.convert_thrust_to_pwm(thrust)

        msg = Float64MultiArray()
        msg.data = [pwm_steer, pwm_thrust, 0.0, 0.0]
        self.publisher_.publish(msg)


    def publish_utm_coordinates(self):
        msg = Float64MultiArray()
        msg.data = [self.x0, self.y0, self.a, self.b]
        self.utm_publisher_.publish(msg)

def main(args=None):
    rclpy.init(args=args)
    node = ActuatorPublisher()
    rclpy.spin(node)
    node.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
