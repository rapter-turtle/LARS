clc; clear; close all;

%% Simulation parameters
dt = 0.01; % Simulation step (fine-grained)
control_update_dt = 0.1; % Control update interval
control_update_steps = control_update_dt / dt; % Control update every 10 steps
T = 16; % Simulation duration
t = 0:dt:T;
N = length(t);

%% Initialize states for x and y
x = zeros(1, N); u = zeros(1, N);
y = zeros(1, N); v = zeros(1, N);
psi = zeros(1, N); r = zeros(1, N);
psi(1) = 70*3.141592/180;
y(1) = -20;
x(1) = 2;

tau_u = zeros(1, N);
tau_r = zeros(1, N);
tau_v = zeros(1, N);

tu = 0;
tr = 0;
tv = 0;

hp = 3.5;
lf = 3.5;
lb = -3;

V0 = 3;
xd_dot = 1*sin(t); % Docking station velocity (x)
yd_dot = V0 + 1*cos(0.5*t); % Docking station velocity (y)
xd = zeros(1, N); % Docking station x position
yd = zeros(1, N); % Docking station y position

x_desired = zeros(1, N);
y_desired = zeros(1, N);
xd_desired = 0;
yd_desired = 0;
ux = 0;
uy = 0;

acceptance_rad = 0.2;

M = 37.758;
I = 18.35; 
Xu = 8.9149;
Xuu = 11.2101;
Nr = 16.9542;
Nrrr = 12.8966;
Yv = 15;
Yvv = 3;
Yr = 6;
Nv = 6;

%% Control gains
Kp = 3;  % Proportional gain
alpha = 0.05;

%% Logging for additional plots
ux_log = zeros(1, N); % Control input log (x)
uy_log = zeros(1, N); % Control input log (y)
safety_log = zeros(1, N); % Safety constraint log


%% Figure Setup
fig = figure('Position', [100, 100, 1200, 600]); % Adjust figure size

% Left: Ship Motion Animation
subplot(3, 4, [1 2 5 6 9 10]); % 왼쪽 큰 플롯 (애니메이션)
hold on;
axis equal;
grid on;
xlabel('X Position');
ylabel('Y Position');
title('Ship Motion Simulation');
% xlim([-20, 20]); ylim([-30, 200]);
xlim([-10, 10]);  % 초기 x 범위는 -30 ~ 5로 설정
ylim([-20, 70]);


% Initialize speed-dependent movement
x_offset = V0 * t;

% Initial ship shape (Pentagon)
ship_x = [3.5 2.5 -3 -3 2.5]; % Shape of the ship (pentagon)
ship_y = [0 1.2 1.2 -1.2 -1.2];

func_x = linspace(-15, 15, 100);
func_plot = plot(func_x, -(func_x - xd(1)).^2 + yd(1), 'k--', 'LineWidth', 1.5);
% Plot initial ship
ship_plot = fill(ship_x, ship_y, 'b'); % Ship shape

% ** Acceptence Radius 업데이트 (선두와 선미 위치에 원 그리기) **
theta = psi(1); % 초기 선박의 회전각
radii = acceptance_rad * ones(1, 2);  % 두 원의 반지름 (선두, 선미)
center_x = [x(1) + lf * cos(theta), x(1) + lb * cos(theta)]; % 선두, 선미의 x 위치
center_y = [y(1) + lf * sin(theta), y(1) + lb * sin(theta)]; % 선두, 선미의 y 위치

% 선두와 선미에 빨간 원 추가
circle_hf = rectangle('Position', [center_x(1)-acceptance_rad, center_y(1)-acceptance_rad, 2*acceptance_rad, 2*acceptance_rad], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
circle_hb = rectangle('Position', [center_x(2)-acceptance_rad, center_y(2)-acceptance_rad, 2*acceptance_rad, 2*acceptance_rad], 'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
circle_dock = rectangle('Position', [0, 0, 2*acceptance_rad, 2*acceptance_rad], 'Curvature', [1, 1], 'EdgeColor', 'g', 'LineWidth', 2);


% Add simulation time display (Top-left corner)
time_text = text(-45, 45, 'Time: 0.00 s', 'FontSize', 12, 'Color', 'r');

% Right: 3x2 Subplots
ax(1) = subplot(3,4,3); hold on; title('X Position');
ax(2) = subplot(3,4,4); hold on; title('Y Position');
ax(3) = subplot(3,4,7); hold on; title('Velocity (u, v)');
ax(4) = subplot(3,4,8); hold on; title('Angular Velocity (r)');
ax(5) = subplot(3,4,11); hold on; title('Control Inputs (\tau_u, \tau_r)');
ax(6) = subplot(3,4,12); hold on; title('Heading Angle (\psi)');

% 선을 미리 생성 (plot 업데이트 최적화)
line_x = plot(ax(1), t(1), x(1), 'b');
line_y = plot(ax(2), t(1), y(1), 'r');
line_u = plot(ax(3), t(1), u(1), 'g');
line_v = plot(ax(3), t(1), v(1), 'm');
line_r = plot(ax(4), t(1), r(1), 'k');
line_tau_u = plot(ax(5), t(1), tau_u(1), 'b');
line_tau_r = plot(ax(5), t(1), tau_r(1), 'r');
line_tau_v = plot(ax(5), t(1), tau_v(1), 'g');
line_psi = plot(ax(6), t(1), psi(1), 'k');


% 세로선(Vertical Line) 추가
for i = 1:6
    vline_x(i) = xline(ax(i), t(1), '--k', 'LineWidth', 1.2);
end

%% Video writer setup
video = VideoWriter('ship_simulation.mp4', 'MPEG-4'); % MP4 format
video.FrameRate = 70; % Frame rate
open(video); % Open video writer to start saving

%% Simulation loop
for i = 2:N
    if mod(i, control_update_steps) == 0 
        %% System variable
        u_dot_cbf = (-Xu*u(i-1) - Xuu*sqrt(u(i-1)*u(i-1))*u(i-1))/M;
        v_dot_cbf = (-Yv*v(i-1) -Yr*r(i-1) - Yvv*sqrt(v(i-1)*v(i-1))*v(i-1))/M;
        r_dot_cbf = (-Nv*v(i-1) -Nr*r(i-1) - Nrrr*r(i-1)*r(i-1)*r(i-1))/I;

        eta_x_dotdot = u_dot_cbf*cos(psi(i-1)) - v_dot_cbf*sin(psi(i-1)) - r_dot_cbf*hp*sin(psi(i-1)) - u(i-1)*r(i-1)*sin(i-1) - v(i-1)*r(i-1)*cos(i-1) - r(i-1)*r(i-1)*hp*cos(i-1);
        eta_y_dotdot = u_dot_cbf*sin(psi(i-1)) + v_dot_cbf*cos(psi(i-1)) + r_dot_cbf*hp*cos(psi(i-1)) + u(i-1)*r(i-1)*cos(i-1) - v(i-1)*r(i-1)*sin(i-1) - r(i-1)*r(i-1)*hp*sin(i-1);
    
        eta_x = x(i-1) + hp*cos(psi(i-1));
        eta_y = y(i-1) + hp*sin(psi(i-1));    
        eta_x_dot = u(i-1)*cos(psi(i-1)) - v(i-1)*sin(psi(i-1)) - r(i-1)*hp*sin(psi(i-1));
        eta_y_dot = u(i-1)*sin(psi(i-1)) + v(i-1)*cos(psi(i-1)) + r(i-1)*hp*cos(psi(i-1));


        %% Define decision variables (X)
        Vd_x = sdpvar(1,1);      
        Vd_dot_x = sdpvar(1,1);
        s = sdpvar(1,1);

        %% Define decision variables (Y)
        Vd_y = sdpvar(1,1);      
        Vd_dot_y = sdpvar(1,1);

        %% Define Cost Function
        Objective = Vd_dot_x^2 + Vd_dot_y^2 + 1e12 * s^2;

        %% Define Constraints (CBF)
        Constraints = [];
        Constraints = [Constraints, -(eta_y - yd(i-1)) - (eta_x - xd(i-1))^2 + alpha * (-(Vd_y - yd_dot(i-1))-2*(eta_x - xd(i-1))*(Vd_x - xd_dot(i-1))) + s >= 0];
        Constraints = [Constraints, s>=0];

        %% Define Constraints (CLF X)
        Constraints = [Constraints, Vd_x == x_desired(i-1) + Vd_dot_x * control_update_dt];
        Constraints = [Constraints, (eta_x - xd(i-1)) * (Vd_x - xd_dot(i-1)) <= 0];
        Constraints = [Constraints, Vd_dot_x <= 2];
        Constraints = [Constraints, Vd_dot_x >= -2];

        %% Define Constraints (CLF Y)
        Constraints = [Constraints, Vd_y == y_desired(i-1) + Vd_dot_y * control_update_dt];
        Constraints = [Constraints, (eta_y - yd(i-1)) * (Vd_y - yd_dot(i-1)) <= 0];
        Constraints = [Constraints, Vd_dot_y <= 2];
        Constraints = [Constraints, Vd_dot_y >= -2];

        %% Solve using SDPT3
        options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        sol = optimize(Constraints, Objective, options);

        xd_desired = value(Vd_dot_x);
        yd_desired = value(Vd_dot_y);
        
        % Compute control input
        ux = -Kp * (eta_x_dot - x_desired(i-1)) + xd_desired - eta_x_dotdot;
        uy = -Kp * (eta_y_dot - y_desired(i-1)) + yd_desired - eta_y_dotdot;

        tu = M*(ux*cos(psi(i-1)) + uy*sin(psi(i-1)));
        tr = I*(uy*cos(psi(i-1)) - ux*sin(psi(i-1)))/hp;

    end
    
    
    %% System update
    % Update states using Euler integration
    xd(i) = xd(i-1) + xd_dot(i) * dt;
    yd(i) = yd(i-1) + yd_dot(i) * dt;

    x_desired(i) = x_desired(i-1) + xd_desired * dt;
    y_desired(i) = y_desired(i-1) + yd_desired * dt;
    
    % Log additional values
    ux_log(i) = ux;
    uy_log(i) = uy;
    safety_log(i) = -(y(i-1) + hp*sin(i-1) - yd(i-1)) - (x(i-1) + hp*cos(i-1) - xd(i-1))^2;


    %% System update

    
    % Update states using Euler integration
    tau_u(i) = tu;
    tau_r(i) = tr;
    tau_v(i) = 0;
    
    x_dot = u(i-1)*cos(psi(i-1)) - v(i-1)*sin(psi(i-1));
    y_dot = u(i-1)*sin(psi(i-1)) + v(i-1)*cos(psi(i-1));
    u_dot = (-Xu*u(i-1) - Xuu*sqrt(u(i-1)*u(i-1))*u(i-1) + tau_u(i))/M;
    v_dot = (-Yv*v(i-1) -Yr*r(i-1) - Yvv*sqrt(v(i-1)*v(i-1))*v(i-1) + tau_v(i))/M;
    r_dot = (-Nv*v(i-1) -Nr*r(i-1) - Nrrr*r(i-1)*r(i-1)*r(i-1) + tau_r(i))/I;
    
    x(i) = x(i-1) + x_dot * dt;
    y(i) = y(i-1) + y_dot * dt;
    psi(i) = psi(i-1) + r(i-1) * dt;
    u(i) = u(i-1) + u_dot * dt;
    v(i) = v(i-1) + v_dot * dt;
    r(i) = r(i-1) + r_dot * dt;

    % ** 실시간 선박 업데이트 **
    theta = psi(i); % 선박의 회전각
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)]; % 회전행렬 적용
    ship_rotated = R * [ship_x; ship_y]; % 회전된 선박 좌표
    
    set(ship_plot, 'XData', ship_rotated(1,:) + x(i), 'YData', ship_rotated(2,:) + y(i));
    % ** 실시간 y 범위 업데이트 (ylim 변경) **
    % y 범위는 V0와 t에 따라 업데이트됨
    current_yrange_min = -30 + V0 * t(i);
    current_yrange_max = 5 + V0 * t(i);
    axis([-10, 10, current_yrange_min, current_yrange_max]);
    drawnow;

    % Update Simulation Time Text
    set(time_text, 'String', sprintf('Time: %.2f s', t(i)));

    % ** Acceptence Radius 업데이트 (선두와 선미 위치에 원 그리기) **
    center_x = [x(i) + lf * cos(psi(i)), x(i) + lb * cos(psi(i))]; % 선두, 선미의 x 위치
    center_y = [y(i) + lf * sin(psi(i)), y(i) + lb * sin(psi(i))]; % 선두, 선미의 y 위치
    
    % 선두와 선미에 빨간 원 추가
    set(circle_hf, 'Position', [center_x(1)-acceptance_rad, center_y(1)-acceptance_rad, 2*acceptance_rad, 2*acceptance_rad]);
    set(circle_hb, 'Position', [center_x(2)-acceptance_rad, center_y(2)-acceptance_rad, 2*acceptance_rad, 2*acceptance_rad]);
    set(circle_dock, 'Position', [xd(i), yd(i), 2*acceptance_rad, 2*acceptance_rad]);

    % Subplot 업데이트 (선을 이어서 그림)
    set(line_x, 'XData', t(1:i), 'YData', x(1:i));
    set(line_y, 'XData', t(1:i), 'YData', y(1:i));
    set(line_u, 'XData', t(1:i), 'YData', u(1:i));
    set(line_v, 'XData', t(1:i), 'YData', v(1:i));
    set(line_r, 'XData', t(1:i), 'YData', r(1:i));
    set(line_tau_u, 'XData', t(1:i), 'YData', tau_u(1:i));
    set(line_tau_r, 'XData', t(1:i), 'YData', tau_r(1:i));
    set(line_tau_v, 'XData', t(1:i), 'YData', tau_v(1:i));
    set(line_psi, 'XData', t(1:i), 'YData', psi(1:i));
    set(func_plot, 'YData', -(func_x - xd(i)).^2 + yd(i));

    % Update vertical lines
    for j = 1:6
        set(vline_x(j), 'Value', t(i));
    end
    % Capture the current frame
    frame = getframe(fig); % Capture the frame from the figure
    writeVideo(video, frame);
    drawnow;
    
    % Stop simulation if figure is closed
    if ~isvalid(fig)
        break;
    end
end
close(video);