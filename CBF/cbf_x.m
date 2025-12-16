clc; clear; close all;
yalmip('clear');

%% Simulation parameters
dt = 0.01; % Simulation step (fine-grained)
control_update_dt = 0.1; % Control update interval
control_update_steps = control_update_dt / dt; % Control update every 10 steps
T = 50; % Simulation duration
t = 0:dt:T;
N = length(t);

%% Initialize states
x = zeros(1, N); x_dot = zeros(1, N);
x(1) = 3;
xd_dot = 1*sin(0.5*t); % Docking station velocity
xd = zeros(1, N); % Docking station x position

x_desired = zeros(1, N);
xd_desired = 0;
ux = 0;

%% Control gains
Kp = 5;  % Proportional gain
alpha = 5;

%% Logging for additional plots
ux_log = zeros(1, N); % Control input log
safety_log = zeros(1, N); % Safety constraint log

%% Simulation loop
for i = 2:N
    % Control input update every 10 steps
    if mod(i, control_update_steps) == 0 
        %% Define decision variables
        Vd = sdpvar(1,1);      % Desired velocity
        Vd_dot = sdpvar(1,1);
        s = sdpvar(1,1);       % Slack variable for constraint relaxation

        %% Define Cost Function
        Objective = Vd_dot^2 + 1000000000*s^2;

        %% Define Constraints
        Constraints = [];
        Constraints = [Constraints, Vd == x_desired(i-1) + Vd_dot * control_update_dt]; % Velocity update equation
        Constraints = [Constraints, (x(i-1) - xd(i-1)) + alpha * (Vd - xd_dot(i-1)) - alpha*(0.3/Kp) + s >= 0]; % Safety constraint
        Constraints = [Constraints, (x(i-1) - xd(i-1))* (Vd - xd_dot(i-1)) <= 0];
        Constraints = [Constraints, Vd_dot <= 1];
        Constraints = [Constraints, Vd_dot >= -1];

        %% Solve using SDPT3
        options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        sol = optimize(Constraints, Objective, options);

        xd_desired = value(Vd_dot);
        
        % Compute control input
        ux = -Kp * (x_dot(i-1) - x_desired(i-1)) + xd_desired;
    end
    
    % Update states using Euler integration
    x_dot(i) = x_dot(i-1) + ux * dt + 0.3*sin(t(i-1))*dt;
    x(i) = x(i-1) + x_dot(i) * dt;
    xd(i) = xd(i-1) + xd_dot(i) * dt;
    x_desired(i) = x_desired(i-1) + xd_desired * dt;
    
    % Log additional values
    ux_log(i) = ux;
    safety_log(i) = (x(i-1) - xd(i-1)) + alpha * (x_desired(i-1) - xd_dot(i-1));
end

%% Plot results
figure;
subplot(2,2,1); % Main trajectory plot
plot(t, xd, 'r--', 'LineWidth', 2); hold on;
plot(t, x, 'b', 'LineWidth', 2);
legend('Docking Station', 'Robot');
xlabel('Time (s)'); ylabel('X Position (m)');
title('Robot Docking Simulation');
grid on;

subplot(2,2,2); % Control input plot
plot(t, ux_log, 'g', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Control Input u_x');
title('Control Input Over Time');
grid on;

subplot(2,2,3); % Desired vs actual docking position plot
plot(t, x_dot, 'r--', 'LineWidth', 2); hold on;
plot(t, x_desired, 'b', 'LineWidth', 2);
legend('Docking Station (x_d)', 'Desired Position (x_{desired})');
xlabel('Time (s)'); ylabel('X Position (m)');
title('Docking Station vs. Desired Position');
grid on;

subplot(2,2,4); % Safety constraint plot
plot(t, safety_log, 'm', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Safety Constraint Value');
title('Safety Constraint Over Time');
grid on;
