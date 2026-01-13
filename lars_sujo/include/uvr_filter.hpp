// uvr_filter.hpp
#pragma once

#include <array>
#include <cmath>
#include <algorithm>

namespace uvr_filter
{

// ============================================
// 공통 유틸
// ============================================
inline double wrapToPi(double a)
{
    while (a >  M_PI) a -= 2.0 * M_PI;
    while (a < -M_PI) a += 2.0 * M_PI;
    return a;
}

inline double unwrapYaw(double yaw_new, double yaw_prev)
{
    double d = yaw_new - yaw_prev;
    if (d >  M_PI) yaw_new -= 2.0 * M_PI;
    if (d < -M_PI) yaw_new += 2.0 * M_PI;
    return yaw_new;
}

// ============================================
// 1) 간단 필터: (x,y,psi) -> (u,v,r)
//    - finite difference로 vx,vy 계산
//    - body frame로 회전하여 u,v 계산
//    - r은 yaw 미분
//    - LPF 적용
// ============================================
struct UVRFilterLPF
{
    bool init = false;

    double x_prev = 0.0;
    double y_prev = 0.0;
    double psi_prev = 0.0;

    // filtered outputs
    double u_hat = 0.0;
    double v_hat = 0.0;
    double r_hat = 0.0;

    // tuning: time constants [s]
    double tau_u = 0.1;
    double tau_v = 0.1;
    double tau_r = 0.1;

    void reset(double x, double y, double psi)
    {
        init = true;
        x_prev = x;
        y_prev = y;
        psi_prev = psi;
        u_hat = v_hat = r_hat = 0.0;
    }

    static inline double alpha(double dt, double tau)
    {
        tau = std::max(1e-3, tau);
        return std::exp(-dt / tau);
    }

    // 입력: global pose (x,y,psi), dt
    // 출력: u_hat, v_hat, r_hat
    void update(double x, double y, double psi_meas, double dt)
    {
        if (dt <= 1e-4) return;

        if (!init)
        {
            reset(x, y, psi_meas);
            return;
        }

        double psi = unwrapYaw(psi_meas, psi_prev);

        // global velocity
        const double vx = (x - x_prev) / dt;
        const double vy = (y - y_prev) / dt;

        // body velocity
        const double c = std::cos(psi);
        const double s = std::sin(psi);

        const double u =  c * vx + s * vy;
        const double v = -s * vx + c * vy;

        // yaw rate
        const double dpsi = wrapToPi(psi - psi_prev);
        const double r = dpsi / dt;

        // LPF
        const double au = alpha(dt, tau_u);
        const double av = alpha(dt, tau_v);
        const double ar = alpha(dt, tau_r);

        u_hat = au * u_hat + (1.0 - au) * u;
        v_hat = av * v_hat + (1.0 - av) * v;
        r_hat = ar * r_hat + (1.0 - ar) * r;

        // store
        x_prev = x;
        y_prev = y;
        psi_prev = psi;
    }
};

// ============================================
// 2) 운동학 EKF: 상태 [x,y,psi,u,v,r]
//    측정 [x,y,psi]만으로 u,v,r 추정
// ============================================
struct UVREKF
{
    bool init = false;

    // state: [x, y, psi, u, v, r]
    std::array<double, 6> x{};

    // covariance 6x6, row-major
    std::array<double, 36> P{};

    // ---- tuning parameters ----
    // measurement noise (분산)
    double R_x   = 0.25;                         // (m^2)  -> std 0.5m
    double R_y   = 0.25;                         // (m^2)
    double R_psi = (2.0 * M_PI / 180.0) * (2.0 * M_PI / 180.0); // (rad^2) -> 2deg

    // process noise (random walk strength) : 분산/초 느낌으로 dt 곱해서 반영
    double Q_u = 4.0;                           // (m/s)^2
    double Q_v = 4.0;                           // (m/s)^2
    // double Q_r = (20.0 * M_PI / 180.0) * (5.0 * M_PI / 180.0);   // (rad/s)^2
    double Q_r = std::pow(20.0 * M_PI/180.0, 2);  // 기존 5deg -> 20deg
    
    void setIdentityP(double diag)
    {
        std::fill(P.begin(), P.end(), 0.0);
        for (int i = 0; i < 6; i++) P[i * 6 + i] = diag;
    }

    void reset(double x_m, double y_m, double psi_m)
    {
        init = true;
        x = {x_m, y_m, psi_m, 0.0, 0.0, 0.0};
        setIdentityP(10.0);
    }

    static inline double& Pij(std::array<double,36>& P, int i, int j) { return P[i * 6 + j]; }
    static inline double  Pij(const std::array<double,36>& P, int i, int j) { return P[i * 6 + j]; }

    // 예측 단계
    void predict(double dt)
    {
        if (!init || dt <= 1e-4) return;

        const double X   = x[0];
        const double Y   = x[1];
        const double psi = x[2];
        const double u   = x[3];
        const double v   = x[4];
        const double r   = x[5];

        const double c = std::cos(psi);
        const double s = std::sin(psi);

        // state propagation
        x[0] = X + dt * (u * c - v * s);
        x[1] = Y + dt * (u * s + v * c);
        x[2] = wrapToPi(psi + dt * r);
        // u,v,r : random-walk (값 자체는 그대로, 공분산에 Q만 추가)

        // Jacobian F (6x6)
        double F[6][6] = {};
        for (int i = 0; i < 6; i++) F[i][i] = 1.0;

        F[0][2] = dt * (-u * s - v * c);
        F[0][3] = dt * ( c);
        F[0][4] = dt * (-s);

        F[1][2] = dt * ( u * c - v * s);
        F[1][3] = dt * ( s);
        F[1][4] = dt * ( c);

        F[2][5] = dt;

        // P = F P F^T + Q
        std::array<double,36> FP{};
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < 6; k++) sum += F[i][k] * Pij(P, k, j);
                FP[i * 6 + j] = sum;
            }
        }

        std::array<double,36> Pnew{};
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < 6; k++) sum += FP[i * 6 + k] * F[j][k]; // *F^T
                Pnew[i * 6 + j] = sum;
            }
        }
        P = Pnew;

        // process noise (dt 스케일)
        Pij(P, 3, 3) += Q_u * dt;
        Pij(P, 4, 4) += Q_v * dt;
        Pij(P, 5, 5) += Q_r * dt;
    }

    // 업데이트 단계 (측정: x,y,psi)
    void update(double x_m, double y_m, double psi_m)
    {
        if (!init)
        {
            reset(x_m, y_m, psi_m);
            return;
        }

        // psi unwrap around current estimate
        psi_m = unwrapYaw(psi_m, x[2]);

        // innovation y = z - h(x)
        double ytilde[3] = {
            x_m - x[0],
            y_m - x[1],
            wrapToPi(psi_m - x[2])
        };

        // S = HPH^T + R (3x3) -> 상단 3x3
        double S[3][3] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                S[i][j] = Pij(P, i, j);

        S[0][0] += R_x;
        S[1][1] += R_y;
        S[2][2] += R_psi;

        // S^{-1} (3x3) via Gauss-Jordan
        double A[3][6] = {
            {S[0][0], S[0][1], S[0][2], 1,0,0},
            {S[1][0], S[1][1], S[1][2], 0,1,0},
            {S[2][0], S[2][1], S[2][2], 0,0,1}
        };

        for (int col = 0; col < 3; col++)
        {
            double piv = A[col][col];
            if (std::fabs(piv) < 1e-12) return;

            for (int j = 0; j < 6; j++) A[col][j] /= piv;

            for (int row = 0; row < 3; row++)
            {
                if (row == col) continue;
                double f = A[row][col];
                for (int j = 0; j < 6; j++) A[row][j] -= f * A[col][j];
            }
        }

        double Sinv[3][3] = {
            {A[0][3], A[0][4], A[0][5]},
            {A[1][3], A[1][4], A[1][5]},
            {A[2][3], A[2][4], A[2][5]}
        };

        // K = P H^T S^{-1} = P(:,0:2) * Sinv  => 6x3
        double K[6][3] = {};
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double sum = 0.0;
                for (int k = 0; k < 3; k++)
                    sum += Pij(P, i, k) * Sinv[k][j];
                K[i][j] = sum;
            }
        }

        // state update: x = x + K*y
        for (int i = 0; i < 6; i++)
        {
            x[i] += K[i][0] * ytilde[0] + K[i][1] * ytilde[1] + K[i][2] * ytilde[2];
        }
        x[2] = wrapToPi(x[2]);

        // covariance update: P = P - K*(H*P)
        // H*P 는 P(0:2,:)
        double HP[3][6] = {};
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 6; j++)
                HP[i][j] = Pij(P, i, j);

        std::array<double,36> Pnew = P;
        for (int i = 0; i < 6; i++)
        {
            for (int j = 0; j < 6; j++)
            {
                double sub = K[i][0] * HP[0][j] + K[i][1] * HP[1][j] + K[i][2] * HP[2][j];
                Pnew[i * 6 + j] = Pij(P, i, j) - sub;
            }
        }
        P = Pnew;
    }

    // 편의 함수: pose measurement 한 번으로 예측+업데이트까지
    void step(double x_m, double y_m, double psi_m, double dt)
    {
        if (!init) { reset(x_m, y_m, psi_m); return; }
        predict(dt);
        update(x_m, y_m, psi_m);
    }

    // 추정된 u,v,r 반환
    double u() const { return x[3]; }
    double v() const { return x[4]; }
    double r() const { return x[5]; }
};

} // namespace uvr_filter
