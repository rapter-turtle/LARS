import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64MultiArray
from sensor_msgs.msg import Imu
import numpy as np
import pyproj


class ShipSimulator(Node):
    def __init__(self):
        super().__init__('ship_simulator')
        
        # ROS2 I/O
        self.state_pub = self.create_publisher(Float64MultiArray, '/ekf/estimated_state', 10)
        self.lars_state_pub = self.create_publisher(Float64MultiArray, '/lars_state', 10)
        self.imu_pub = self.create_publisher(Imu, '/imu/data', 10)
        self.control_sub = self.create_subscription(Float64MultiArray, '/actuator_outputs', self.control_callback, 10)

        # 초기 상태 [x, y, ψ, u, v, r, δ, F_eff]
        self.ship_state = np.zeros(6)
        self.ship_state[0] = (289577.66 + 291591.05) * 0.5 - 25.0# X
        self.ship_state[1] = (4117065.30 + 4118523.52) * 0.5 - 25.0 # Y
        self.ship_state[2] = 0 * np.pi / 180.0  # ψ
        self.disturbance_state = np.zeros((6))

        # 시뮬레이션 step
        self.dt = 0.01  
        self.con_dt = 0.5  

        self.now_delta = 0.0
        self.now_F = 0.0

        self.nn = 0
        self.tt = 0

        # 입력 지연 설정
        self.delay = 3.0   # 2초 delay
        self.delay_step = int(self.delay / self.con_dt)
        self.input_buffer = np.zeros((self.delay_step, 2))  # [δ_cmd, F_cmd] 히스토리
        self.delayed_input = np.array([0,0])

        # Timer
        self.timer = self.create_timer(self.dt, self.run)

        # 좌표 변환
        self.utm_proj = pyproj.Proj(proj="utm", zone=52, datum="WGS84")
        self.wgs84_proj = pyproj.Proj(proj="latlong", datum="WGS84")

        self.MS_x = (289577.66 + 291591.05) * 0.5
        self.MS_y = (4117065.30 + 4118523.52) * 0.5
        self.MS_psi = 0 * np.pi / 180.0

        self.CD_x = (289577.66 + 291591.05) * 0.5
        self.CD_y = (4117065.30 + 4118523.52) * 0.5 - 10.0
        self.CD_psi = 0 * np.pi / 180.0

        self.stream_speed = 3.0


    def wave_disturbance(self, disturbance_state, wave_direction, wind_speed, omega, lamda, Kw, sigmaF1, sigmaF2, dt):
        omega_e = np.abs(omega - (omega * omega / 9.81) * wind_speed * np.cos(wave_direction))
        x1 = disturbance_state[0]
        x2 = disturbance_state[1]

        omegaF1 = np.random.normal(0.0, sigmaF1)
        omegaF2 = np.random.normal(0.0, sigmaF2)

        xdot = np.array([x2, -omega_e * omega_e * x1 - 2 * lamda * omega_e * x2 + Kw * omegaF1, omegaF2])
        disturbance_state = xdot * dt + disturbance_state

        disturbance_force = disturbance_state[1] + disturbance_state[2]
        return disturbance_state, disturbance_force
    

    def control_callback(self, msg):
        """actuator_outputs → [steer_pwm, thrust_pwm, bow, 0]"""
        steer_pwm = msg.data[0]
        thrust_pwm = msg.data[1]

        # PWM → δ, F 변환
        delta_cmd = self.convert_pwm_to_steering(steer_pwm)
        F_cmd = self.convert_pwm_to_thrust(thrust_pwm)
        self.now_F = F_cmd
        self.now_delta = delta_cmd
        self.delayed_input = [delta_cmd, F_cmd]


    def convert_pwm_to_steering(self, pwm):
        # print(pwm)
        if pwm >= 2000.0: return 300.0
        elif pwm >= 1550.0 and pwm < 2000.0: return (pwm - 1550.0) / 1.6667
        elif pwm <= 1450.0 and pwm > 1000.0: return (pwm - 1450.0) / 1.6667
        elif pwm <= 1000.0: return -300.0
        else: return 0.0


    def convert_pwm_to_thrust(self, pwm):
        """Inverse of convert_thrust_to_pwm: PWM → rpm_thrust (with deadzone)."""

        # Deadzone center
        if 1450.0 < pwm < 1550.0:
            return 0.0

        # Negative thrust region
        elif pwm <= 1450.0:
            thrust = (pwm - 1450.0) / 3.9   # linear thrust
            thrust_2 = -(thrust**2)         # 복원된 thrust_2 (음수 영역)
            rpm_thrust = 0.0#thrust_2 - (22.0**2)  # deadzone 보정

        # Positive thrust region
        elif pwm >= 1550.0:
            thrust = (pwm - 1550.0) / 3.9   # linear thrust
            thrust_2 = (thrust**2)          # 복원된 thrust_2 (양수 영역)
            rpm_thrust = thrust_2  # deadzone 보정

        else:
            rpm_thrust = 0.0

        return 0.01*rpm_thrust

    def utm_to_latlon(self, utm_x, utm_y):
        """Convert UTM coordinates to latitude and longitude."""
        lon, lat = pyproj.transform(self.utm_proj, self.wgs84_proj, utm_x, utm_y)
        return lat, lon

    def ship_dynamics(self, state, control, dt):
        """컨트롤러 모델과 맞춘 선박 dynamics"""

        # Run Parameters
        M, I = 1.0, 1.0

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

        eps = 1e-6

        # States
        x, y, psi, u, v, r = state
        delta = self.now_delta
        F_cmd = self.now_F
        

        # # ------------------------------------- Disturbance ------------------------------------- 
        wind_direction = 0.1
        wind_speed = 3.0

        self.disturbance_state[:3], XY_wave_force = self.wave_disturbance(
            self.disturbance_state[:3], wind_direction, wind_speed, 0.8, 0.1, 0.64, 6, 2, dt)
        self.disturbance_state[3:6], N_wave_force = self.wave_disturbance(
            self.disturbance_state[3:6], wind_direction, wind_speed, 0.8, 0.1, 1.0, 1, 0.1, dt)

        X_wave_force = XY_wave_force * np.cos(wind_direction - psi)
        Y_wave_force = XY_wave_force * np.sin(wind_direction - psi)

        U_wave_force = X_wave_force * np.cos(psi) + Y_wave_force * np.sin(psi)
        V_wave_force = -X_wave_force * np.sin(psi) + Y_wave_force * np.cos(psi)

        # Wind disturbance
        Afw = 0.3
        Alw = 0.44
        LOA = 1.3
        lau = 1.2
        CD_lAF = 0.55
        delta1 = 0.6
        CDl = CD_lAF * Afw / Alw
        CDt = 0.8

        u_rel_wind = u - wind_speed * np.cos(wind_direction - psi)
        v_rel_wind = v - wind_speed * np.sin(wind_direction - psi)
        gamma = -np.arctan2(v_rel_wind, u_rel_wind)

        Cx = CD_lAF * np.cos(gamma) / (1 - delta1 * 0.5 * (1 - CDl / CDt) * (np.sin(2 * gamma))**2)
        Cy = CDt * np.sin(gamma) / (1 - delta1 * 0.5 * (1 - CDl / CDt) * (np.sin(2 * gamma))**2)
        Cn = -0.18 * (gamma) * Cy
        # Cn = -0.18 * (gamma - np.pi * 0.5) * Cy

        X_wind_force = 0.02 * lau * (u_rel_wind**2 + v_rel_wind**2) * Cx * Afw
        Y_wind_force = 0.02 * lau * (u_rel_wind**2 + v_rel_wind**2) * Cy * Alw
        N_wind_force = 0.005 * lau * (u_rel_wind**2 + v_rel_wind**2) * Cn * Alw * LOA

        U_wind_force = X_wind_force * np.cos(psi) + Y_wind_force * np.sin(psi)
        V_wind_force = -X_wind_force * np.sin(psi) + Y_wind_force * np.cos(psi)

        u_dis = (U_wind_force + U_wave_force)*0.4
        v_dis = (V_wind_force + V_wave_force)*0.4
        r_dis = N_wind_force
        
        if u_dis > 0.4:
            u_dis = 0.4
        elif u_dis < -0.4:
            u_dis = -0.4

        if v_dis > 0.4:
            v_dis = 0.4
        elif v_dis < -0.4:
            v_dis = -0.4

        
        # # ------------------------------------- Disturbance -------------------------------------
        

        
        # Dynamics
        u_dot = (- Xu*u - Xuu * np.sqrt(u * u + eps) * u + alpha1*F_cmd*np.cos(alpha2*delta))/M + u_dis
        v_dot = (-Yv*v - Yr*r - Yvv * np.sqrt(v* v + eps) * v + alpha3*F_cmd*np.sin(alpha2*delta))/M + v_dis
        r_dot = (- Nr*r - Nv*v - Nrr * np.sqrt(r * r + eps) * r - alpha4*F_cmd*np.sin(alpha2*delta))/I + r_dis
        # print(u, v, r, F_cmd)

        # Kinematics
        x_dot = u*np.cos(psi) - v*np.sin(psi) - self.stream_speed
        y_dot = u*np.sin(psi) + v*np.cos(psi)
        psi_dot = r

        # Forward Euler integration
        next_state = np.zeros_like(state)
        next_state[0] = x + dt*x_dot
        next_state[1] = y + dt*y_dot
        next_state[2] = psi + dt*psi_dot
        next_state[3] = u + dt*u_dot
        next_state[4] = v + dt*v_dot
        next_state[5] = r + dt*r_dot


        # wrap heading
        if next_state[2] > np.pi: next_state[2] -= 2*np.pi
        if next_state[2] < -np.pi: next_state[2] += 2*np.pi

        # ----- CBF -----
        # # Side CBF
        # d = 4.0
        # x_prime = self.MS_x + d*np.cos(self.MS_psi + 0.5*3.141592)
        # y_prime = self.MS_y + d*np.sin(self.MS_psi + 0.5*3.141592)
        # h_side = y - np.tan(self.MS_psi)*x - y_prime + np.tan(self.MS_psi)*x_prime 
        # # print(h_side)

        # # Dock CBF
        # hp = 3.0
        # a = 3.0
        # b = 0.5 
        # x0 = -5.0
        # y0 = 4.0
        # xh = x + hp*np.cos(psi)    
        # yh = y + hp*np.sin(psi)
        # xr = xh - self.CD_x    
        # yr = yh - self.CD_y
        # xh_dot = -self.stream_speed + u*np.cos(psi) - v*np.sin(psi) - hp*r*np.sin(psi)
        # yh_dot = u*np.sin(psi) + v*np.cos(psi) + hp*r*np.cos(psi)

        # h_dock_left = a*np.tanh(b*(x0 - (np.cos(self.CD_psi)*xr + np.sin(self.CD_psi)*yr))) + y0 - np.sin(self.CD_psi)*xr + np.cos(self.CD_psi)*yr
        # h_dock_right = a*np.tanh(b*(x0 - (np.cos(self.CD_psi)*xr + np.sin(self.CD_psi)*yr))) + y0 + np.sin(self.CD_psi)*xr - np.cos(self.CD_psi)*yr
        # # print("[Dock CBF] Left : ",h_dock_left,", ",h_dock_right )


        return next_state

    def run(self):
        """메인 시뮬레이션 루프"""
  
        self.ship_state = self.ship_dynamics(self.ship_state, self.delayed_input, self.dt)
        lat, lon = self.utm_to_latlon(self.ship_state[0], self.ship_state[1])

        # Prepare the state message to publish
        state_msg = Float64MultiArray()
        # Now send lat, lon along with the ship's state (psi, u, v, r)
        state_msg.data = np.concatenate((self.ship_state,self.ship_state)).tolist()
        # state_msg.data = np.concatenate(([lat, lon], self.ship_state[2:],[0,0,0,0,0,0,0,0,0,0,0])).tolist()

        lars_state_msg = Float64MultiArray()
        lars_state_msg.data = [self.CD_x, self.CD_y, self.CD_psi, self.MS_x, self.MS_y, self.MS_psi, self.stream_speed]
        self.lars_state_pub.publish(lars_state_msg)


        if self.nn >= 10:
            self.state_pub.publish(state_msg)
            self.nn = 0

        self.nn += 1
        self.tt += 1

        

        self.publish_imu()



    def publish_imu(self):
        imu_msg = Imu()
        imu_msg.header.stamp = self.get_clock().now().to_msg()
        imu_msg.header.frame_id = "base_link"
        imu_msg.angular_velocity.z = self.ship_state[5]
        imu_msg.linear_acceleration.x = self.ship_state[3]
        imu_msg.linear_acceleration.y = self.ship_state[4]
        self.imu_pub.publish(imu_msg)

def main(args=None):
    rclpy.init(args=args)
    simulator = ShipSimulator()
    rclpy.spin(simulator)
    simulator.destroy_node()
    rclpy.shutdown()

if __name__ == '__main__':
    main()
