clc; clear; close all;
yalmip('clear');

%% Given values (known parameters)
Vd_0 = 1.0;   % Initial desired velocity
dt = 0.1;     % Time step
x0 = 5.0;     % Initial position
xt = 2.0;     % Target position
alpha = 1.5;  % Safety coefficient
xt_dot = 0.5; % Target velocity

%% Define decision variables
Vd = sdpvar(1,1);      % Desired velocity
Vd_dot = sdpvar(1,1);  % Acceleration (control variable)

%% Define Cost Function
Objective = Vd_dot^2;

%% Define Constraints
Constraints = [];
Constraints = [Constraints, Vd == Vd_0 + Vd_dot * dt]; % Velocity update equation
Constraints = [Constraints, (x0 - xt) + alpha * (Vd - xt_dot) >= 0]; % Safety constraint

%% Solve using SDPT3
options = sdpsettings('solver', 'sdpt3', 'verbose', 1);
sol = optimize(Constraints, Objective, options);

%% Display Results
if sol.problem == 0
    disp('Optimal Solution:');
    disp(['Vd = ', num2str(value(Vd))]);
    disp(['Vd_dot = ', num2str(value(Vd_dot))]);
    disp(['Optimal Cost: ', num2str(value(Objective))]);
else
    disp('QP Solver failed.');
    disp(sol.info);
end
