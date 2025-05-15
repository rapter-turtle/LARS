clear; clc;

% 데이터 불러오기
data = readtable('C:\Users\user\Desktop\ToolboxLS\Examples\Reachability\output_data_Nx71_Kp1_Kd1.csv');

% 필터링 조건: VFlat < 0, xFlat < -5
filtered = data(data.VFlat < 0.0 & data.xFlat < -2, :);

% 조건에 맞는 샘플이 100개 이상 있는지 확인
num_samples = height(filtered);
num_trials = min(10000, num_samples);  % 최대 100개

if num_trials < 1
    error('조건에 맞는 데이터가 충분하지 않습니다.');
end

% 100개 무작위 추출
rng(0);  % 재현 가능성을 위한 시드 고정
rand_indices = randperm(num_samples, num_trials);
init_samples = filtered(rand_indices, :);

% 시뮬레이션 파라미터
dt = 0.01;
T = 6.5;
N = T / dt;
V0 = 2.5;

% 시각화 초기화
figure; hold on;

% 🔹조건에 맞는 점들 먼저 표시 (하늘색 점)
all_vflat_neg = data(data.VFlat < 0.0, :);
scatter(all_vflat_neg.xFlat, all_vflat_neg.yFlat, 40, 'c', 'filled', ...
    'DisplayName', 'All VFlat < 0 Points');


plot(0, 0, 'bo', 'MarkerSize', 10, 'DisplayName', 'Target (xt, yt)');

% 회색 직사각형 추가
% 아래쪽 박스
rectangle('Position', [-1, -2, 2, 0.5], 'FaceColor', [0.2 0.2 0.2 1.0]);

% 위쪽 박스
rectangle('Position', [-1, 1.5, 2, 0.5], 'FaceColor', [0.2 0.2 0.2 1.0]);

% 오른쪽 박스
rectangle('Position', [1, -2, 0.5, 4], 'FaceColor', [0.2 0.2 0.2 1.0]);


rectangle('Position', [-1, -1.5, 2, 3], ...
    'EdgeColor', 'k', 'FaceColor', [0.8 0.8 0.8 0.5]);

% 모든 초기 조건에 대해 시뮬레이션 실행
for i = 1:num_trials
    % 초기 상대 상태 추출
    init_x = init_samples.xFlat(i);
    init_y = init_samples.yFlat(i);
    init_z = init_samples.zFlat(i);  % xo_dot
    init_r = init_samples.rFlat(i);  % yo_dot

    % 절대 위치로 변환 (target 시작점은 [5, 0])
    xo0 = 5 + init_x;
    yo0 = 0 + init_y;
    xo_dot0 = init_z;
    yo_dot0 = init_r;
    state = zeros(6, N+1);
    state(:,1) = [5; 0; xo0; yo0; xo_dot0; yo_dot0];

    for k = 1:N
        % 외란 생성
        Vx = -1 + 2*rand();
        Vy = -1 + 2*rand();
        wx = -0.3 + 0.6*rand();
        wy = -0.3 + 0.6*rand();

        % 현재 상태
        xt = state(1,k); yt = state(2,k);
        xo = state(3,k); yo = state(4,k);
        xo_dot = state(5,k); yo_dot = state(6,k);

        % 속도 업데이트
        xt_dot = V0 + Vx;
        yt_dot = Vy;

        % 제어 입력
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

        state(:,k+1) = [xt_next; yt_next; xo_next; yo_next; xo_dot_next; yo_dot_next];
    end


    % 상대 위치 계산
    rel_x = state(3,:) - state(1,:);
    rel_y = state(4,:) - state(2,:);
    
    % x > 0인 부분 제거
    valid_idx = rel_x <= -0.5;
    rel_x = rel_x(valid_idx);
    rel_y = rel_y(valid_idx);
    
    % 필터링 후 플롯
    plot(rel_x, rel_y, 'r', 'LineWidth', 1, 'HandleVisibility', 'off');

end

% 플롯 설정
axis equal; grid on;
xlabel('x_r [m]'); ylabel('y_r [m]');
xlim([-10,1.5]);
% legend;
% title('Object Trajectories Relative to Moving Target (100 Trials)');

% 제어 함수
function [ux, uy] = control_input(xo, yo, xo_dot, yo_dot, xt_dot, yt_dot, xt, yt)
    kp = 1.0;
    kd = 1.0;
    ux = kp * (xt - xo) + kd * (xt_dot - xo_dot);
    uy = kp * (yt - yo) + kd * (yt_dot - yo_dot);
end
