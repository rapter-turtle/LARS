% clear; clc;

% 데이터 불러오기
% data = readtable('C:\Users\user\Desktop\ToolboxLS\Examples\Reachability\reachability_plot\output_data_Nx71_Kp1_Kd1.csv');

% 필터링 조건: VFlat < 0, xFlat < -5
filtered = data(data.VFlat < 0.0 & data.xFlat < -2, :);

% 조건에 맞는 샘플이 100개 이상 있는지 확인
num_samples = height(filtered);
num_trials = min(1000, num_samples);  % 최대 100개

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
hp = 1;

% 시각화 초기화
figure; hold on;


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
    init_x = init_samples.xFlat(i) - hp;
    init_y = init_samples.yFlat(i);
    init_z = init_samples.zFlat(i);  % xo_dot
    init_r = init_samples.rFlat(i);  % yo_dot

    % 절대 위치로 변환 (target 시작점은 [5, 0])
    xo0 = 5 + init_x;
    yo0 = 0 + init_y;
    xo_dot0 = init_z;
    yo_dot0 = init_r;
    state = zeros(8, N+1);
    hp_state = zeros(2, N+1);
    hp_state(:,1) = [xo0;yo0];
    state(:,1) = [5; 0; xo0; yo0; xo_dot0; yo_dot0; 0; 0];


    for k = 1:N
        % 외란 생성
        Vx = -1 + 2*rand();
        Vy = -1 + 2*rand();
        wx = -0.2 + 0.4*rand();
        wy = -0.2 + 0.4*rand();
        wr = -0.1 + 0.2*rand();

        % 현재 상태
        xt = state(1,k); yt = state(2,k);
        xo = state(3,k); yo = state(4,k);
        uo = state(5,k); vo = state(6,k);
        psio = state(7,k); ro = state(8,k);

        % 속도 업데이트
        xt_dot = V0 + Vx;
        yt_dot = Vy;

        x_hp = xo + hp*cos(psio);
        y_hp = yo + hp*sin(psio);

        x_dot_hp = uo*cos(psio) - vo*sin(psio) - ro*hp*sin(psio);
        y_dot_hp = uo*sin(psio) + vo*cos(psio) + ro*hp*cos(psio);

        % 제어 입력
        [uu, ur] = control_input(x_hp, y_hp, x_dot_hp, y_dot_hp, xt_dot, yt_dot, xt, yt, uo, vo, psio, ro);


        uo_dot = (-uo*Xu - Xuu*sqrt(uo*uo)*uo + uu)/M + wx;
        vo_dot = (-vo*Yv - Yvv*sqrt(vo*vo)*vo - Yr*ro)/M + wy;
        ro_dot = (-vo*Nv - Nr*ro - Nrrr*ro*ro*ro + ur)/I + wr;


        xo_dot = uo*cos(psio) - vo*sin(psio);
        yo_dot = uo*sin(psio) + vo*cos(psio);


        % 적분
        xt_next = xt + xt_dot * dt;
        yt_next = yt + yt_dot * dt;
        uo_next = uo + uo_dot * dt;
        vo_next = vo + vo_dot * dt;
        ro_next = ro + ro_dot * dt;
        xo_next = xo + xo_dot * dt;
        yo_next = yo + yo_dot * dt;
        psio_next = psio + ro_next * dt;

        if psio_next > pi
            psio_next = psio_next - 2*pi;
        elseif psio_next <-pi
            psio_next = psio_next + 2*pi;
        end

        state(:,k+1) = [xt_next; yt_next; xo_next; yo_next; uo_next; vo_next; psio_next; ro_next];
    end


    % 상대 위치 계산
    xhp = state(3,:) + hp*cos(state(7,:));
    yhp = state(4,:) + hp*sin(state(7,:));
    rel_x = xhp - state(1,:);
    rel_y = yhp - state(2,:);
    
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
function [uu, ur] = control_input(xo, yo, xo_dot, yo_dot, xt_dot, yt_dot, xt, yt, u, v, psi, r)
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
    hp = 1;

    uo_dot = (-u*Xu - Xuu*sqrt(u*u)*u)/M;
    vo_dot = (-v*Yv - Yvv*sqrt(v*v)*v - Yr*r)/M;
    ro_dot = (-v*Nv - Nr*r - Nrrr*r*r*r)/I;

    kp = 1.0;
    kd = 1.0;
    ux = kp * (xt - xo) + kd * (xt_dot - xo_dot) - (uo_dot*cos(psi) - vo_dot*sin(psi) - ro_dot*hp*sin(psi));
    uy = kp * (yt - yo) + kd * (yt_dot - yo_dot) - (uo_dot*sin(psi) + vo_dot*cos(psi) + ro_dot*hp*cos(psi));


    uu = (ux*cos(psi) + uy*sin(psi))*M;
    ur = (uy*cos(psi) - ux*sin(psi))*I/hp;

end
