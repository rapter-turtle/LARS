clc; clear; close all;
yalmip('clear');

%% Simulation parameters
dt = 0.01; % Simulation step (fine-grained)
control_update_dt = 0.1; % Control update interval
control_update_steps = control_update_dt / dt; % Control update every 10 steps
T = 40; % Simulation duration
t = 0:dt:T;
N = length(t);

%% Initialize states for x and y
x = zeros(1, N); x_dot = zeros(1, N);
y = zeros(1, N); y_dot = zeros(1, N);
y(1) = -10; 

xd_dot = 1*sin(0.5*t); % Docking station velocity (x)
yd_dot = 1*cos(0.5*t); % Docking station velocity (y)
xd = zeros(1, N); % Docking station x position
yd = zeros(1, N); % Docking station y position

x_desired = zeros(1, N);
y_desired = zeros(1, N);
xd_desired = 0;
yd_desired = 0;
ux = 0;
uy = 0;

%% Control gains
Kp = 5;  % Proportional gain
alpha = 5;

%% Logging for additional plots
ux_log = zeros(1, N); % Control input log (x)
uy_log = zeros(1, N); % Control input log (y)
safety_log = zeros(1, N); % Safety constraint log

%% Simulation loop
for i = 2:N
    if mod(i, control_update_steps) == 0 
        %% Define decision variables (X)
        Vd_x = sdpvar(1,1);      
        Vd_dot_x = sdpvar(1,1);
        s = sdpvar(1,1);

        %% Define decision variables (Y)
        Vd_y = sdpvar(1,1);      
        Vd_dot_y = sdpvar(1,1);

        %% Define Cost Function
        Objective = Vd_dot_x^2 + Vd_dot_y^2 + 1e9 * s^2;

        %% Define Constraints (CBF)
        Constraints = [];
        Constraints = [Constraints, -(y(i-1) - yd(i-1)) - (x(i-1) - xd(i-1))^2 + alpha * (-(Vd_y - yd_dot(i-1))-2*(x(i-1) - xd(i-1))*(Vd_x - xd_dot(i-1))) - alpha*(0.3/Kp) - 2*alpha*abs((x(i-1) - xd(i-1)))*(0.3/Kp) + s >= 0];
        Constraints = [Constraints, s>=0];

        %% Define Constraints (X)
        Constraints = [Constraints, Vd_x == x_desired(i-1) + Vd_dot_x * control_update_dt];
        Constraints = [Constraints, (x(i-1) - xd(i-1)) * (Vd_x - xd_dot(i-1)) <= 0];
        Constraints = [Constraints, Vd_dot_x <= 1];
        Constraints = [Constraints, Vd_dot_x >= -1];

        %% Define Constraints (Y)
        Constraints = [Constraints, Vd_y == y_desired(i-1) + Vd_dot_y * control_update_dt];
        Constraints = [Constraints, (y(i-1) - yd(i-1)) * (Vd_y - yd_dot(i-1)) <= 0];
        Constraints = [Constraints, Vd_dot_y <= 1];
        Constraints = [Constraints, Vd_dot_y >= -1];

        %% Solve using SDPT3
        options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        sol = optimize(Constraints, Objective, options);

        xd_desired = value(Vd_dot_x);
        yd_desired = value(Vd_dot_y);
        
        % Compute control input
        ux = -Kp * (x_dot(i-1) - x_desired(i-1)) + xd_desired;
        uy = -Kp * (y_dot(i-1) - y_desired(i-1)) + yd_desired;
    end
    
    % Update states using Euler integration
    x_dot(i) = x_dot(i-1) + ux * dt + 0.3*sin(t(i-1))*dt;
    y_dot(i) = y_dot(i-1) + uy * dt + 0.3*cos(t(i-1))*dt;
    
    x(i) = x(i-1) + x_dot(i) * dt;
    y(i) = y(i-1) + y_dot(i) * dt;

    xd(i) = xd(i-1) + xd_dot(i) * dt;
    yd(i) = yd(i-1) + yd_dot(i) * dt;

    x_desired(i) = x_desired(i-1) + xd_desired * dt;
    y_desired(i) = y_desired(i-1) + yd_desired * dt;
    
    % Log additional values
    ux_log(i) = ux;
    uy_log(i) = uy;
    safety_log(i) = -(y(i-1) - yd(i-1)) - (x(i-1) - xd(i-1))^2;
end

%% Plot results (Original 6 Subplots)
figure;
subplot(2,3,1);
plot(xd, yd, 'r--', 'LineWidth', 2); hold on;
plot(x, y, 'b', 'LineWidth', 2);
legend('Docking Station', 'Robot');
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Robot Docking Simulation (XY Trajectory)');
grid on;

subplot(2,3,2);
plot(t, ux_log, 'g', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Control Input u_x');
title('Control Input Over Time (X)');
grid on;

subplot(2,3,3);
plot(t, uy_log, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Control Input u_y');
title('Control Input Over Time (Y)');
grid on;

subplot(2,3,4);
plot(t, x_dot, 'r--', 'LineWidth', 2); hold on;
plot(t, x_desired, 'b', 'LineWidth', 2);
legend('Docking Station (x_d)', 'Desired Position (x_{desired})');
xlabel('Time (s)'); ylabel('X Position (m)');
title('Docking Station vs. Desired Position (X)');
grid on;

subplot(2,3,5);
plot(t, y_dot, 'r--', 'LineWidth', 2); hold on;
plot(t, y_desired, 'b', 'LineWidth', 2);
legend('Docking Station (y_d)', 'Desired Position (y_{desired})');
xlabel('Time (s)'); ylabel('Y Position (m)');
title('Docking Station vs. Desired Position (Y)');
grid on;

subplot(2,3,6);
plot(t, safety_log, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Safety Constraint Value');
title('Safety Constraint Over Time');
grid on;

%% Animation Setup
video_filename = 'docking_simulation.mp4';
video_writer = VideoWriter(video_filename, 'MPEG-4');
video_writer.FrameRate = 30;
open(video_writer);

figure;
hold on;
grid on;
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Robot Docking Animation');
xlim([-15 15]); ylim([-15 15]);
plot_handle = plot(x(1), y(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
xd_yd_handle = plot(xd(1), yd(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
func_x = linspace(-15, 15, 100);
func_plot = plot(func_x, -(func_x - xd(1)).^2 + yd(1), 'k--', 'LineWidth', 1.5);

for i = 1:10:N
    set(plot_handle, 'XData', x(i), 'YData', y(i));
    set(xd_yd_handle, 'XData', xd(i), 'YData', yd(i));
    set(func_plot, 'YData', -(func_x - xd(i)).^2 + yd(i));
    frame = getframe(gcf);
    writeVideo(video_writer, frame);
    pause(0.01);
end

close(video_writer);
disp(['Animation saved as ', video_filename]);
