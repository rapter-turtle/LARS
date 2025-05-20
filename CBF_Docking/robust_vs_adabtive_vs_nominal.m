clc; clear; close all;
yalmip('clear');

%% Simulation parameters
dt = 0.01; % Simulation step (fine-grained)
control_update_dt = 0.1; % Control update interval
control_update_steps = control_update_dt / dt; % Control update every 10 steps
T = 100; % Simulation duration
t = 0:dt:T;
N = length(t);

%% Initialize states for x and y
x = zeros(1, N); x_dot = zeros(1, N);
y = zeros(1, N); y_dot = zeros(1, N);
y(1) = -10; 

xd_dot = 0*sin(0.5*t); % Docking station velocity (x)
yd_dot = 0*cos(0.5*t); % Docking station velocity (y)
xd = 0*ones(1, N); % Docking station x position
yd = -1.0*ones(1, N); % Docking station y position

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
DOB_bound = 0.1;

w_max = 0.2;
% wx = w_max*sin(1.0*t);
% wy = w_max*cos(1.0*t);
wx = 0.5*w_max*ones(1, N);
wy = 0.5*w_max*ones(1, N);

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
Kp = 0.05;  % Proportional gain
Kd = 2;
alpha1 = 10;
alpha2 = 10;
alpha_l = 0.1;
alpha_r = 0.1;
a = 2;

x(1) = -2;
y(1) = -10;
%% Logging for additional plots
ux_log = zeros(1, N); % Control input log (x)
uy_log = zeros(1, N); % Control input log (y)
safety_log1 = zeros(1, N); % Safety constraint log
safety_log2 = zeros(1, N); % Safety constraint log


%% Setup figure for animation + 3x2 subplot layout (with improved spacing)
hFig = figure('Name', 'Robot Docking + Live Subplots', 'Position', [100, 100, 1600, 900]);

% Left: Large XY Trajectory
ax1 = axes('Position', [0.05, 0.1, 0.5, 0.8]);
hDock = plot(ax1, xd(1), yd(1), 'ro', 'MarkerSize', 5); hold(ax1, 'on');
hRobot = plot(ax1, x(1), y(1), 'bo', 'MarkerSize', 5);
x_line = linspace(-15, 15, 100);
plot(ax1, x_line, a*x_line, 'k--', 'LineWidth', 1.5);
plot(ax1, x_line, -a*x_line, 'k--', 'LineWidth', 1.5);
legend(ax1, {'Docking Station', 'Robot'}, 'Location', 'best');
xlabel(ax1, 'X Position (m)'); ylabel(ax1, 'Y Position (m)');
title(ax1, 'Robot Docking Simulation (XY Trajectory)');
axis(ax1, 'equal'); xlim(ax1, [-15 15]); ylim(ax1, [-15 2]); grid(ax1, 'on');

% Define right-side layout: 3x2 grid with more vertical spacing
right_left = 0.6;
right_width = 0.35;
col_gap = 0.04;
row_gap = 0.07;
plot_width = (right_width - col_gap) / 2;
plot_height = (0.8 - 2 * row_gap) / 3;
base_y_top = 0.87;  % slightly raised

% Row 1
ax2 = axes('Position', [right_left, base_y_top - plot_height, plot_width, plot_height]);
hUx = animatedline(ax2, 'Color', 'g', 'LineWidth', 2);
xlabel(ax2, 'Time (s)'); ylabel(ax2, 'u_x'); title(ax2, 'Control Input X'); grid(ax2, 'on');

ax3 = axes('Position', [right_left + plot_width + col_gap, base_y_top - plot_height, plot_width, plot_height]);
hUy = animatedline(ax3, 'Color', 'b', 'LineWidth', 2);
xlabel(ax3, 'Time (s)'); ylabel(ax3, 'u_y'); title(ax3, 'Control Input Y'); grid(ax3, 'on');

% Row 2
ax4 = axes('Position', [right_left, base_y_top - 1*(plot_height + row_gap) - plot_height, plot_width, plot_height]);
hCBF1 = animatedline(ax4, 'Color', 'b', 'LineWidth', 2);
xlabel(ax4, 'Time (s)'); ylabel(ax4, 'Left CBF'); title(ax4, 'Left CBF'); grid(ax4, 'on');

ax5 = axes('Position', [right_left + plot_width + col_gap, base_y_top - 1*(plot_height + row_gap) - plot_height, plot_width, plot_height]);
hCBF2 = animatedline(ax5, 'Color', 'b', 'LineWidth', 2);
xlabel(ax5, 'Time (s)'); ylabel(ax5, 'Right CBF'); title(ax5, 'Right CBF'); grid(ax5, 'on');

% Row 3
ax6 = axes('Position', [right_left, base_y_top - 2*(plot_height + row_gap) - plot_height, plot_width, plot_height]);
hWx = animatedline(ax6, 'Color', 'r', 'LineWidth', 2);
hWxEst = animatedline(ax6, 'Color', 'b', 'LineWidth', 2);
xlabel(ax6, 'Time (s)'); ylabel(ax6, 'x disturbance'); title(ax6, 'Disturbance'); ylim(ax6, [-w_max w_max]); grid(ax6, 'on');

ax7 = axes('Position', [right_left + plot_width + col_gap, base_y_top - 2*(plot_height + row_gap) - plot_height, plot_width, plot_height]);
% hReserved = plot(ax7, 0, 0, 'k:');
hWy = animatedline(ax7, 'Color', 'r', 'LineWidth', 2);
hWyEst = animatedline(ax7, 'Color', 'b', 'LineWidth', 2);
xlabel(ax7, 'Time (s)'); ylabel(ax7, 'Reserved');
title(ax7, 'Custom Slot (TODO)'); grid(ax7, 'on');






%% Simulation loop
for i = 2:N
    if mod(i, control_update_steps) == 0 
  
        %% Original control input

        % ux_origin = 0.0;%-Kp * (x(i-1) - xd(i-1)) - Kd * (x_dot(i-1) - xd_dot(i-1));
        % uy_origin = 0.0;%-Kp * (y(i-1) - yd(i-1)) - Kd * (y_dot(i-1) - yd_dot(i-1));

        ux_origin = -Kp * (x(i-1) - xd(i-1)) - Kd * (x_dot(i-1) - xd_dot(i-1));
        uy_origin = -Kp * (y(i-1) - yd(i-1)) - Kd * (y_dot(i-1) - yd_dot(i-1));

        % ux_origin = -Kp * (x(i-1) - xd(i-1));
        % uy_origin = 2;

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
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (a*ux_opt - uy_opt) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*ux_opt - uy_opt) + s2 >= 0];
       
        %% Nominal Robust
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (a*ux_opt - uy_opt) - (a*w_max + w_max) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*ux_opt - uy_opt) - (a*w_max + w_max) + s2 >= 0];
              
        %% fea_Nominal
        Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s1 >= 0];
        Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max)) + s2 >= 0];
         
        %% fea_RaCBF #1
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) + (a*x_dis_dot - y_dis_dot) - (a*DOB_bound + DOB_bound) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*DOB_bound + DOB_bound) + (a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) + (-a*x_dis_dot - y_dis_dot) - (a*DOB_bound + DOB_bound) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*DOB_bound + DOB_bound) + (-a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s2 >= 0];
        
        %% fea_RaCBF #2
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) + (a*x_dis_dot - y_dis_dot) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*DOB_bound + DOB_bound) + (a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) + (-a*x_dis_dot - y_dis_dot) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*DOB_bound + DOB_bound) + (-a*x_dis_estim(i-1) - y_dis_estim(i-1))) + s2 >= 0];
                
        %% fea_RCBF
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) -(a*w_max + w_max) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*w_max + w_max)) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) -(a*w_max + w_max) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (-a*w_max + w_max)) + s2 >= 0];

        %% fixed fea_RCBF
        % Constraints = [Constraints, (alpha1 + alpha2)*(a*ux_opt - uy_opt) + alpha1*alpha2*(a*x_dot(i-1) - y_dot(i-1)) -(a*w_max + w_max) + alpha_l*((alpha1 + alpha2)*(a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (a*w_max + w_max)) + s1 >= 0];
        % Constraints = [Constraints, (alpha1 + alpha2)*(-a*ux_opt - uy_opt) + alpha1*alpha2*(-a*x_dot(i-1) - y_dot(i-1)) -(a*w_max + w_max) + alpha_r*((alpha1 + alpha2)*(-a*x_dot(i-1) - y_dot(i-1)) + alpha1*alpha2*(-a*x(i-1) - y(i-1)) + (-a*f_max - f_max) - (-a*w_max + w_max)) + s2 >= 0];
        % Constraints = [Constraints, uy_opt >= a*ux_opt + (-a*f_max - f_max) + (a*w_max + w_max)];
        % Constraints = [Constraints, uy_opt >= -a*ux_opt + (-a*f_max - f_max) + (a*w_max + w_max)];
        % Constraints = [Constraints, uy_opt <= a*ux_opt - (-a*f_max - f_max) - (a*w_max + w_max)];
        % Constraints = [Constraints, uy_opt <= -a*ux_opt - (-a*f_max - f_max) - (a*w_max + w_max)];

        %% Define Constraints (ux, uy)
        Constraints = [Constraints, ux_opt <= f_max];
        Constraints = [Constraints, ux_opt >= -f_max];
        Constraints = [Constraints, uy_opt <= f_max];
        Constraints = [Constraints, uy_opt >= -f_max];
        % 
        Constraints = [Constraints, ux_opt <= ux_log(i-1) + 0.5*control_update_dt*f_max];
        Constraints = [Constraints, ux_opt >= ux_log(i-1) - 0.5*control_update_dt*f_max];
        Constraints = [Constraints, uy_opt <= uy_log(i-1) + 0.5*control_update_dt*f_max];
        Constraints = [Constraints, uy_opt >= uy_log(i-1) - 0.5*control_update_dt*f_max];

        %% Solve using SDPT3
        options = sdpsettings('solver', 'sdpt3', 'verbose', 0);
        sol = optimize(Constraints, Objective, options);

        ux = value(ux_opt);
        uy = value(uy_opt);
        value(s1)
        value(s2)
  

    end
    
    % Update states using Euler integration
    x_dot(i) = x_dot(i-1) + ux * dt + wx(i-1) * dt;
    y_dot(i) = y_dot(i-1) + uy * dt + wy(i-1) * dt;
    
    x(i) = x(i-1) + x_dot(i) * dt;
    y(i) = y(i-1) + y_dot(i) * dt;


    if mod(i, control_update_steps) == 0 
        %% DOB
        x_error_k = x_estim_k - x_dot(i);
        y_error_k = y_estim_k - y_dot(i); 
        x_estim_k = x_estim(i-1) + control_update_dt*(ux + x_dis_estim(i-1) + As*x_error_k);
        y_estim_k = y_estim(i-1) + control_update_dt*(uy + y_dis_estim(i-1) + As*y_error_k);
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
    if mod(i,10) == 0
        set(hDock, 'XData', xd(i), 'YData', yd(i));
        set(hRobot, 'XData', x(i), 'YData', y(i));

        addpoints(hUx, t(i), ux_log(i));
        addpoints(hUy, t(i), uy_log(i));
        addpoints(hCBF1, t(i), safety_log1(i));
        addpoints(hCBF2, t(i), safety_log2(i));
        addpoints(hWx, t(i), wx(i));
        addpoints(hWxEst, t(i), x_dis_estim(i));
        addpoints(hWy, t(i), wy(i));
        addpoints(hWyEst, t(i), y_dis_estim(i));

        drawnow;
        pause(0.01);
    end

end

% %% Plot results (Original 6 Subplots)
% figure;
% subplot(2,3,1);
% plot(xd, yd, 'r--', 'LineWidth', 2); hold on;
% plot(x, y, 'b', 'LineWidth', 2);
% legend('Docking Station', 'Robot');
% xlabel('X Position (m)'); ylabel('Y Position (m)');
% title('Robot Docking Simulation (XY Trajectory)');
% grid on;
% 
% subplot(2,3,2);
% plot(t, ux_log, 'g', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Control Input u_x');
% title('Control Input Over Time (X)');
% grid on;
% 
% subplot(2,3,3);
% plot(t, uy_log, 'b', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('Control Input u_y');
% title('Control Input Over Time (Y)');
% grid on;
% 
% subplot(2,3,4);
% plot(t, safety_log1, 'b', 'LineWidth', 2); hold on;
% % legend('Docking Station (x_d)', 'Desired Position (x_{desired})');
% xlabel('Time (s)'); ylabel('Left CBF');
% title('Left CBF');
% grid on;
% 
% subplot(2,3,5);
% plot(t, safety_log2, 'b', 'LineWidth', 2);
% % legend('Docking Station (y_d)', 'Desired Position (y_{desired})');
% xlabel('Time (s)'); ylabel('Right CBF');
% title('Right CBF');
% grid on;
% 
% subplot(2,3,6);
% plot(t, wx, 'r', 'LineWidth', 2); hold on;
% plot(t, x_dis_estim, 'b', 'LineWidth', 2);
% xlabel('Time (s)'); ylabel('x');
% ylim([-w_max w_max]);
% title('Disturbance');
% grid on;
