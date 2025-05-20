clc; clear; close all;
yalmip('clear');

%% Simulation parameters
dt = 0.01; % Simulation step (fine-grained)
control_update_dt = 0.1; % Control update interval
control_update_steps = control_update_dt / dt; % Control update every 10 steps
T = 20; % Simulation duration
t = 0:dt:T;
N = length(t);

%% Initialize states for x and y
x = zeros(1, N); x_dot = zeros(1, N);
y = zeros(1, N); y_dot = zeros(1, N);
y(1) = -10; 

xd_dot = 0*sin(0.5*t); % Docking station velocity (x)
yd_dot = 0*cos(0.5*t); % Docking station velocity (y)
xd = 0*ones(1, N); % Docking station x position
yd = -1*ones(1, N); % Docking station y position

x_estim_k = 0;
y_estim_k = 0;
x_error_k = 0;
y_error_k = 0;
x_dis_estim_k = 0;
y_dis_estim_k = 0;
        
x_error = zeros(1, N); 
y_error = zeros(1, N);
x_estim = zeros(1, N);
y_estim = zeros(1, N);
x_dis_estim = zeros(1, N);
y_dis_estim = zeros(1, N);
x_dis_dot = zeros(1, N);
y_dis_dot = zeros(1, N);
As = -1;
DOB_bound = 0.05;

w_max = 0.3;
wx = w_max*sin(1.0*t);
wy = w_max*cos(1.0*t);

x_desired = zeros(1, N);
y_desired = zeros(1, N);
xd_desired = 0;
yd_desired = 0;
ux = 0;
uy = 0;
ux_origin = 0;
uy_origin = 0;

f_max = 1;

%% Control gains
Kp = 0.5;  % Proportional gain
Kd = 2;
alpha1 = 10;
alpha2 = 10;
alpha_l = 0.5;
alpha_r = 0.5;
a = 1;

x(1) = -5;
y(1) = -10;
%% Logging for additional plots
ux_log = zeros(1, N); % Control input log (x)
uy_log = zeros(1, N); % Control input log (y)
safety_log1 = zeros(1, N); % Safety constraint log
safety_log2 = zeros(1, N); % Safety constraint log



%% Setup figure for animation
figure('Name', 'Robot Docking Animation');
hDock = plot(xd(1), yd(1), 'ro', 'MarkerSize', 5, 'DisplayName', 'Docking Station'); hold on;
hRobot = plot(x(1), y(1), 'bo', 'MarkerSize', 5, 'DisplayName', 'Robot');

% y = a*x, y = -a*x 직선 그리기
x_line = linspace(-15, 15, 100);
y_left = a * x_line;
y_right = -a * x_line;
hLineLeft = plot(x_line, y_left, 'k--', 'LineWidth', 1.5, 'DisplayName', 'y = a x');
hLineRight = plot(x_line, y_right, 'k--', 'LineWidth', 1.5, 'DisplayName', 'y = -a x');

legend('Location', 'best');
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('Robot Docking Simulation (XY Trajectory)');
axis equal;
xlim([-15 15]); ylim([-15 2]);
grid on;



%% Simulation loop
for i = 2:N
    if mod(i, control_update_steps) == 0 
  
        %% Original control input

        ux_origin = 0;%-Kp * (x(i-1) - xd(i-1)) - Kd * (x_dot(i-1) - xd_dot(i-1));
        uy_origin = -Kp * (y(i-1) - yd(i-1)) - Kd * (y_dot(i-1) - yd_dot(i-1));

        %% Define decision variables (X)
        ux_opt = sdpvar(1,1);      
        uy_opt = sdpvar(1,1);
        s1 = sdpvar(1,1);
        s2 = sdpvar(1,1);

        %% Define Cost Function
        Objective = (ux_opt - ux_origin)^2 + (uy_opt - uy_origin)^2 + 1e9 * s1^2 + 1e9 * s2^2;

        %% Define Constraints (CBF)
        Constraints = [];
        %% Nominal
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s2 >= 0];
        
        %% aCBF
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) + (a*x_dis_dot - y_dis_dot) - (a*DOB_bound + DOB_bound) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max) + (a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) + (-a*x_dis_dot - y_dis_dot) - (a*DOB_bound + DOB_bound) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max) + (-a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s2 >= 0];
        
        %% RCBF
        Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) - (a*w_max + w_max) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s1 >= 0];
        Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) - (a*w_max + w_max) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) + y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s2 >= 0];
        

        % Constraints = [Constraints, s1>=0];
        % Constraints = [Constraints, s2>=0];

        %% Define Constraints (ux, uy)
        Constraints = [Constraints, ux_opt <= f_max];
        Constraints = [Constraints, ux_opt >= -f_max];
        Constraints = [Constraints, uy_opt <= f_max];
        Constraints = [Constraints, uy_opt >= -f_max];

        %% Solve using SDPT3
        options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        sol = optimize(Constraints, Objective, options);

        ux = value(ux_opt);
        uy = value(uy_opt);
        % ux = 0;
        % uy = 0;
        

        %% DOB
        x_estim_k = x_estim(i-1) + control_update_dt*(ux + x_dis_estim(i-1) + As*x_error(i-1));
        y_estim_k = y_estim(i-1) + control_update_dt*(uy + y_dis_estim(i-1) + As*y_error(i-1));
        x_error_k = x_estim_k - x_dot(i-1);
        y_error_k = y_estim_k - y_dot(i-1); 
        x_dis_estim_k = -As*exp(As*control_update_dt)*x_error_k/(exp(As*control_update_dt) - 1);
        y_dis_estim_k = -As*exp(As*control_update_dt)*y_error_k/(exp(As*control_update_dt) - 1);
        
        if x_dis_estim_k >= w_max
            x_dis_estim_k = w_max;
        elseif x_dis_estim_k < -w_max
            x_dis_estim_k = -w_max;
        end
        if y_dis_estim_k >= w_max
            y_dis_estim_k = w_max;
        elseif y_dis_estim_k < -w_max
            y_dis_estim_k = -w_max;
        end      

    end
    
    % Update states using Euler integration
    x_dot(i) = x_dot(i-1) + ux * dt + wx(i-1) * dt;
    y_dot(i) = y_dot(i-1) + uy * dt + wy(i-1) * dt;
    
    x(i) = x(i-1) + x_dot(i) * dt;
    y(i) = y(i-1) + y_dot(i) * dt;

    % Log additional values
    ux_log(i) = ux;
    uy_log(i) = uy;
    safety_log1(i) = a*x_dot(i-1) - y_dot(i-1) + alpha1*(a*x(i-1) - y(i-1));
    safety_log2(i) = -a*x_dot(i-1) - y_dot(i-1) + alpha1*(-a*x(i-1) - y(i-1));
    x_error(i) = x_error_k;
    y_error(i) = y_error_k;
    x_estim(i) = x_estim_k;
    y_estim(i) = y_estim_k;
    x_dis_estim(i) = x_dis_estim_k;
    y_dis_estim(i) = y_dis_estim_k;
    x_dis_dot(i) = As*x_error(i);
    y_dis_dot(i) = As*y_error(i);

%% animation
    if mod(i,5) == 0
        set(hDock, 'XData', xd(i), 'YData', yd(i));
        set(hRobot, 'XData', x(i), 'YData', y(i));
        drawnow;
        pause(0.01); % Optional: 조절 가능
    end

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
plot(t, safety_log1, 'b', 'LineWidth', 2); hold on;
% legend('Docking Station (x_d)', 'Desired Position (x_{desired})');
xlabel('Time (s)'); ylabel('Left CBF');
title('Left CBF');
grid on;

subplot(2,3,5);
plot(t, safety_log2, 'b', 'LineWidth', 2);
% legend('Docking Station (y_d)', 'Desired Position (y_{desired})');
xlabel('Time (s)'); ylabel('Right CBF');
title('Right CBF');
grid on;

subplot(2,3,6);
plot(t, wx, 'r', 'LineWidth', 2); hold on;
plot(t, x_dis_estim, 'b', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('x');
ylim([-w_max w_max]);
title('Disturbance');
grid on;

