% simulation.m
clear; clc;

% 시뮬레이션 파라미터
dt = 0.01;         % 시간 간격
T = 6.5;            % 총 시뮬레이션 시간 (초)
N = T / dt;        % 스텝 수
V0 = 2.5;          % 상수 전방 속도

% 상태: [xt, yt, xo, yo, xo_dot, yo_dot]
state = zeros(6, N+1);
state(:,1) = [5; 0; 0; 0; 0; 0];  % 초기 상태 설정

% 기록용
time = 0:dt:T;


% 메인 시뮬레이션 루프
for k = 1:N
    % 외란 생성
    Vx = -1 + 2*rand();           % [-1, 1]
    Vy = -1 + 2*rand();           % [-1, 1]
    wx = -0.3 + 0.6*rand();       % [-0.3, 0.3]
    wy = -0.3 + 0.6*rand();       % [-0.3, 0.3]

    % 현재 상태
    xt = state(1,k);
    yt = state(2,k);
    xo = state(3,k);
    yo = state(4,k);
    xo_dot = state(5,k);
    yo_dot = state(6,k);

    % 상태 업데이트
    xt_dot = V0 + Vx;
    yt_dot = Vy;

    % 제어 입력 (원하면 이 함수 수정)
    [ux, uy] = control_input(xo, yo, xo_dot, yo_dot, xt_dot, yt_dot, xt, yt);

    xo_dot_dot = ux + wx;
    yo_dot_dot = uy + wy;

    % 적분
    xt_next = xt + xt_dot * dt;
    yt_next = yt + yt_dot * dt;
    xo_dot_next = xo_dot + xo_dot_dot * dt;
    yo_dot_next = yo_dot + yo_dot_dot * dt;
    xo_next = xo + xo_dot * dt;
    yo_next = yo + yo_dot * dt;

    % 다음 상태 저장
    state(:,k+1) = [xt_next; yt_next; xo_next; yo_next; xo_dot_next; yo_dot_next];
end

% 상대 좌표 계산
rel_x = state(3,:) - state(1,:);  % xo - xt
rel_y = state(4,:) - state(2,:);  % yo - yt

% 시각화
figure;
plot(0, 0, 'bo', 'MarkerSize', 10, 'DisplayName', 'Target (xt, yt)'); hold on;
plot(rel_x, rel_y, 'r', 'LineWidth', 2, 'DisplayName', 'Object relative to Target');

% 회색 직사각형 추가: x [-1,1], y [-1.5,1.5]
rectangle('Position', [-1, -1.5, 2, 3], 'EdgeColor', 'k', 'FaceColor', [0.8 0.8 0.8 0.5]);

axis equal; grid on;
xlabel('Relative X'); ylabel('Relative Y');
legend;
title('Object Trajectory Relative to Target');


function [ux, uy] = control_input(xo, yo, xo_dot, yo_dot, xt_dot, yt_dot, xt, yt)
    % 단순 추적 제어 (P 제어)
    kp = 1.0;
    kd = 1.0;
    ux = kp * (xt - xo) + kd * (xt_dot - xo_dot);
    uy = kp * (yt - yo) + kd * (yt_dot - yo_dot);
end

