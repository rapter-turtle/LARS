clc; clear;

% 계수 정의
a1 = 5;
a2 = 5;
a_l = 5.0;

% 그리드 정의
x = linspace(-1, 1, 400);
y = linspace(-5, 5, 400);
[X, Y] = meshgrid(x, y);

% 부등식 영역
f0 = - Y;
f1 = -a1.*X - Y;
f2 = -1 - (a1 + a2).*Y - a1*a2.*X;
f3 = 1*(a1 + a2) - a_l - (a1*a2 + a_l*(a1 + a2)).*Y - a_l*a1*a2.*X;

% 조건 만족 영역
region0 = f0 >= 0;
region1 = f1 >= 0;
region2 = f2 >= 0;
region_both = region0 & region1 & region2;

% region_both에 해당하는 (x, y, f3)만 추출
X_region = X(region_both);
Y_region = Y(region_both);
F3_region = f3(region_both);

% 3D plot
figure;
scatter3(X_region, Y_region, F3_region, 10, F3_region, 'filled');
xlabel('x'); ylabel('y'); zlabel('f3(x, y)');
title('f3 in the Feasible Region');
colorbar;
grid on;
view(45, 30);
