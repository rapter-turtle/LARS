#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import csv
import time
import math
import threading
from datetime import datetime

import rclpy
from rclpy.node import Node
from std_msgs.msg import Float64MultiArray
from pynput import keyboard


def clamp(x, lo, hi):
    return max(lo, min(hi, x))


def yaw_unwrap(yaw_new, yaw_prev):
    """
    yaw 래핑(-pi~pi) 값이 들어와도 연속적으로 풀어서(unwrap) 저장.
    """
    if yaw_prev is None:
        return yaw_new
    d = yaw_new - yaw_prev
    if d > math.pi:
        yaw_new -= 2.0 * math.pi
    elif d < -math.pi:
        yaw_new += 2.0 * math.pi
    return yaw_new


class KeyboardControlNode(Node):
    def __init__(self):
        super().__init__('keyboard_control_node_usv')

        # =========================
        # (1) 토픽 설정
        # =========================
        self.actuator_pub = self.create_publisher(Float64MultiArray, '/actuator_outputs', 10)

        self.ekf_sub = self.create_subscription(
            Float64MultiArray,
            '/ekf/estimated_state',
            self.ekf_callback,
            10
        )

        # =========================
        # (2) 조작/제어 파라미터
        # =========================
        self.dt = 0.1  # 타이머 주기 (초)

        # 입력 범위 (당신이 원하는 UI 조작 범위)
        self.STEER_UI_MIN = -300.0
        self.STEER_UI_MAX =  300.0

        # 실제 조향 라디안 범위: ±(pi*60/180) = ±pi/3
        self.STEER_RAD_MAX = math.pi * 60.0 / 180.0  # 1.0472...

        # 추력 범위
        self.THRUST_MIN = 0.0
        self.THRUST_MAX = 5000.0

        # 키 입력 1회당 목표값 변화량(“목표값”에 반영)
        self.STEER_STEP_UI = 10.0     # 좌/우 한번 누를 때 목표 조향(UI) 변화량
        self.THRUST_STEP = 50.0       # 위/아래 한번 누를 때 목표 추력 변화량

        # 조향 자동 복귀(매 타이머마다 목표 조향을 0쪽으로 끌어당김)
        self.STEER_RETURN_PER_TICK = 8.0  # dt=0.1이면 초당 80 UI 정도로 중앙 복귀

        # ===== 변화율 제한 (rate limit) =====
        # 타이머(dt)마다 적용될 최대 변화량
        self.MAX_DSTEER_UI_PER_TICK = 5.0   # dt=0.1 -> 초당 300 UI
        self.MAX_DTHRUST_PER_TICK   = 200.0  # dt=0.1 -> 초당 2000 thrust

        # =========================
        # (3) 상태 변수
        # =========================
        self.lock = threading.Lock()

        # EKF 상태 (msg order: x,y,yaw,u,v,r)
        self.x = 0.0
        self.y = 0.0
        self.yaw = 0.0
        self.yaw_prev = None
        self.u = 0.0
        self.v = 0.0
        self.r = 0.0
        self.have_state = False

        # 사용자 목표 입력(키보드로 바꾸는 값)
        self.target_steer_ui = 0.0   # -300~300
        self.target_thrust   = 0.0   # 0~5000
        self.bow = 0.0              # q/w/e로 -1/0/1

        # 실제 적용 입력(변화율 제한 후 퍼블리시하는 값)
        self.applied_steer_ui = 0.0
        self.applied_thrust   = 0.0

        # =========================
        # (4) CSV 로깅
        # =========================
        now = datetime.now().strftime("%Y%m%d_%H%M%S")
        default_path = os.path.join(os.path.expanduser("~"), f"usv_keyboard_log_{now}.csv")
        self.csv_path = default_path
        self.log_buffer = []  # 종료 시 한 번에 저장

        self.get_logger().info(f"[LOG] CSV 저장 경로: {self.csv_path}")

        # =========================
        # (5) 타이머 시작
        # =========================
        self.timer = self.create_timer(self.dt, self.timer_callback)

        # =========================
        # (6) 키보드 리스너
        # =========================
        self.listener = keyboard.Listener(on_press=self.on_press)
        self.listener.start()

        self.get_logger().info("키보드 제어 시작:")
        self.get_logger().info("  ↑/↓ : thrust 증가/감소")
        self.get_logger().info("  ←/→ : steer 좌/우")
        self.get_logger().info("  q/w/e : bow -1/0/1")
        self.get_logger().info("  r : steer/thrust 리셋")
        self.get_logger().info("  esc : 종료(리스너 stop) + Ctrl+C로 노드 종료 권장")


    # -----------------------------
    # EKF 콜백: 상태 업데이트
    # -----------------------------
    def ekf_callback(self, msg: Float64MultiArray):
        if len(msg.data) < 6:
            return

        with self.lock:
            self.x = float(msg.data[0])
            self.y = float(msg.data[1])

            yaw_raw = float(msg.data[2])
            self.yaw = yaw_unwrap(yaw_raw, self.yaw_prev)
            self.yaw_prev = self.yaw

            self.u = float(msg.data[3])
            self.v = float(msg.data[4])
            self.r = float(msg.data[5])

            self.have_state = True


    # -----------------------------
    # UI steer -> 실제 라디안 변환
    # -----------------------------
    def steer_ui_to_rad(self, steer_ui: float) -> float:
        steer_ui = clamp(steer_ui, self.STEER_UI_MIN, self.STEER_UI_MAX)
        return (steer_ui / self.STEER_UI_MAX) * self.STEER_RAD_MAX


    # -----------------------------
    # 변화율 제한 적용 (현재 -> 목표)
    # -----------------------------
    def rate_limit(self, current, target, max_delta):
        delta = target - current
        delta = clamp(delta, -max_delta, +max_delta)
        return current + delta


    # -----------------------------
    # 타이머 콜백: 퍼블리시 + 로깅
    # -----------------------------
    def timer_callback(self):
        with self.lock:
            # (A) 목표값에 대한 clamp
            self.target_steer_ui = clamp(self.target_steer_ui, self.STEER_UI_MIN, self.STEER_UI_MAX)
            self.target_thrust   = clamp(self.target_thrust,   self.THRUST_MIN,   self.THRUST_MAX)

            # (B) 조향 자동 복귀(목표값을 0쪽으로 이동)
            if self.target_steer_ui > 0.0:
                self.target_steer_ui -= self.STEER_RETURN_PER_TICK
                if self.target_steer_ui < 0.0:
                    self.target_steer_ui = 0.0
            elif self.target_steer_ui < 0.0:
                self.target_steer_ui += self.STEER_RETURN_PER_TICK
                if self.target_steer_ui > 0.0:
                    self.target_steer_ui = 0.0

            # (C) 변화율 제한 적용해서 “실제 적용 입력” 갱신
            self.applied_steer_ui = self.rate_limit(
                self.applied_steer_ui,
                self.target_steer_ui,
                self.MAX_DSTEER_UI_PER_TICK
            )
            self.applied_thrust = self.rate_limit(
                self.applied_thrust,
                self.target_thrust,
                self.MAX_DTHRUST_PER_TICK
            )

            # (D) 실제 퍼블리시 값(라디안/추력)
            steer_rad = self.steer_ui_to_rad(self.applied_steer_ui)
            thrust    = float(self.applied_thrust)

            # actuator publish
            out = Float64MultiArray()
            # 원하는 포맷으로: [steer(rad), thrust(0~5000), bow, 0.0]
            out.data = [steer_rad, thrust, float(self.bow), 0.0]
            self.actuator_pub.publish(out)

            # (E) 로깅: 상태 + 입력을 매 tick 기록
            # 시간은 ROS time nanosec 사용
            t_ns = self.get_clock().now().nanoseconds
            t_sec = t_ns * 1e-9

            # ====== (추가) 매 타임스텝 제어입력 출력 ======
            # 10Hz로 계속 찍힘
            self.get_logger().info(
                f"[10Hz] u_ctrl: steer_ui={self.applied_steer_ui:+7.2f}, "
                f"steer_rad={steer_rad:+6.3f} rad, "
                f"thrust={thrust:7.1f}, bow={self.bow:+.1f}"
            )

            if self.have_state:
                row = [
                    t_sec,
                    self.x, self.y, self.yaw, self.u, self.v, self.r,
                    self.applied_steer_ui, steer_rad,
                    thrust,
                    float(self.bow),
                    self.target_steer_ui,
                    self.target_thrust
                ]
                self.log_buffer.append(row)


    # -----------------------------
    # 키 입력 처리
    # -----------------------------
    def on_press(self, key):
        try:
            with self.lock:
                if key == keyboard.Key.up:
                    self.target_thrust += self.THRUST_STEP
                elif key == keyboard.Key.down:
                    self.target_thrust -= self.THRUST_STEP
                elif key == keyboard.Key.left:
                    self.target_steer_ui -= self.STEER_STEP_UI
                elif key == keyboard.Key.right:
                    self.target_steer_ui += self.STEER_STEP_UI
                elif key == keyboard.Key.esc:
                    # 리스너만 멈추고, 노드 자체는 Ctrl+C로 종료하는 걸 권장
                    self.listener.stop()
                    self.get_logger().info("ESC 입력: 키보드 리스너를 중지했습니다. (노드 종료는 Ctrl+C)")
                elif key == keyboard.KeyCode.from_char('r'):
                    self.reset_values()
                elif key == keyboard.KeyCode.from_char('q'):
                    self.bow = -1.0
                elif key == keyboard.KeyCode.from_char('w'):
                    self.bow = 0.0
                elif key == keyboard.KeyCode.from_char('e'):
                    self.bow = 1.0

                # clamp는 timer에서 다시 하지만, 여기서도 한 번 더 안전하게
                self.target_steer_ui = clamp(self.target_steer_ui, self.STEER_UI_MIN, self.STEER_UI_MAX)
                self.target_thrust   = clamp(self.target_thrust,   self.THRUST_MIN,   self.THRUST_MAX)

        except Exception:
            pass


    def reset_values(self):
        self.target_steer_ui = 0.0
        self.target_thrust = 0.0
        self.applied_steer_ui = 0.0
        self.applied_thrust = 0.0
        self.bow = 0.0
        self.get_logger().info("리셋: steer/thrust/bow 모두 0으로 초기화")


    # -----------------------------
    # 종료 시 CSV 저장
    # -----------------------------
    def save_csv(self):
        header = [
            "t_sec",
            "x", "y", "yaw", "u", "v", "r",
            "applied_steer_ui", "applied_steer_rad",
            "applied_thrust",
            "bow",
            "target_steer_ui",
            "target_thrust"
        ]

        try:
            os.makedirs(os.path.dirname(self.csv_path), exist_ok=True)
            with open(self.csv_path, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(header)
                writer.writerows(self.log_buffer)

            self.get_logger().info(f"[LOG] CSV 저장 완료: {self.csv_path} (rows={len(self.log_buffer)})")
        except Exception as e:
            self.get_logger().error(f"[LOG] CSV 저장 실패: {e}")


def main(args=None):
    rclpy.init(args=args)
    node = KeyboardControlNode()

    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        # 키보드 리스너 정리
        try:
            node.listener.stop()
        except Exception:
            pass

        # CSV 저장
        node.save_csv()

        node.destroy_node()
        rclpy.shutdown()


if __name__ == "__main__":
    main()
