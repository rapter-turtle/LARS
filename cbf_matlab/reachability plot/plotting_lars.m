clear; clc;

% 데이터 불러오기
data = readtable('C:\Users\user\Desktop\ToolboxLS\Examples\Reachability\output_data_Nx71_Kp1_Kd1.csv');

% 절댓값이 0.01보다 작은 VFlat 값 필터링
filtered = data(abs(data.VFlat) < 0.01, :);

% figure 생성 및 크기 지정
figure('Position', [100, 100, 1600, 600]);

%% --- Subplot 1: (x, y, z) ---
subplot(1, 2, 1);
scatter3(filtered.xFlat, filtered.yFlat, filtered.zFlat, 50, filtered.zFlat, 'filled');
xlabel('x'); ylabel('y'); zlabel('z');
title('$x$, $y$, $\dot{x}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('$x$[m]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$y$[m]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('$\dot{x}$[m/s]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on; colorbar; colormap jet;
axis equal; view(45, 25);

% 회색 박스 추가 (x: -1~1, y: -2~2, z: -5~5)
drawOpaqueBox([-1 1.2], [-2 2], [-5.2 5.2]);

%% --- Subplot 2: (x, y, r) ---
subplot(1, 2, 2);
scatter3(filtered.xFlat, filtered.yFlat, filtered.rFlat, 50, filtered.rFlat, 'filled');
xlabel('x'); ylabel('y'); zlabel('r');
title('$x$, $y$, $\dot{y}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('$x$[m]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('$y$[m]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
zlabel('$\dot{y}$[m/s]', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on; colorbar; colormap jet;
axis equal; view(45, 25);

% 회색 박스 추가 (x: -1~1, y: -2~2, r: -5~5)
drawOpaqueBox([-1 1.2], [-2 2], [-5.2 5.2]);

% 조명
lighting phong; camlight headlight;

%% === 함수: 꽉 찬 불투명 회색 박스 그리기 ===
function drawOpaqueBox(xlim, ylim, zlim)
    % 각 꼭짓점 좌표 계산
    [X, Y, Z] = ndgrid([xlim(1), xlim(2)], [ylim(1), ylim(2)], [zlim(1), zlim(2)]);
    vertices = [X(:), Y(:), Z(:)];

    % 육면체 각 면의 인덱스 정의 (6개 면)
    faces = [
        1 2 4 3;  % bottom (z = zmin)
        5 6 8 7;  % top (z = zmax)
        1 2 6 5;  % side x+
        3 4 8 7;  % side x-
        1 3 7 5;  % side y-
        2 4 8 6   % side y+
    ];

    % patch로 불투명 회색 박스 생성 (모서리 검정)
    patch('Vertices', vertices, 'Faces', faces, ...
          'FaceColor', [0.6 0.6 0.6], ...
          'FaceAlpha', 1.0, ...            % 불투명
          'EdgeColor', 'k', ...            % 모서리 검정
          'LineWidth', 1.2);               % 모서리 두께 선택적
end
