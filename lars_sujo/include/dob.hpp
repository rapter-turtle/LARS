#pragma once

#include <array>
#include <cmath>

namespace DOB
{

// ===============================
// DOB State container
// ===============================
struct DOBState
{
    std::array<double,3> state_estim {{0.0, 0.0, 0.0}};     // [û, v̂, r̂]
    std::array<double,3> param_estim {{0.0, 0.0, 0.0}};     // disturbance estimate
    std::array<double,3> param_filtered {{0.0, 0.0, 0.0}};  // filtered disturbance
};

// ===============================
// Disturbance Observer Update
// ===============================
inline void updateDOB(
    const std::array<double,8>& state,   // [x,y,psi,u,v,r,delta,F]
    DOBState& dob,
    double dt, 
    double stream_speed, 
    double Xu, 
    double Xuu,
    double Yv,
    double Yvv, 
    double Yr, 
    double Yrr, 
    double Yrrr,
    double Nr, 
    double Nrr, 
    double Nrrr, 
    double Nv, 
    double u_F, 
    double v_F, 
    double r_F, 
    double Xu_dot, 
    double Yv_dot, 
    double Nr_dot,
    double M, 
    double I
)
{
    // -----------------------------
    // Model parameters
    // -----------------------------
    // constexpr double Xu = 0.0845;
    // constexpr double Xuu = 0.0195;
    // constexpr double Yv = 0.0485;
    // constexpr double Yvv = 0.0988;
    // constexpr double Yr = 0.151;
    // constexpr double Nr = 0.6939;
    // constexpr double Nrr = 0.0;
    // constexpr double Nv = 0.0;

    // constexpr double alpha1 = 0.0452;
    // constexpr double alpha2 = 2.0 / 500.0;
    // constexpr double alpha3 = 0.0188;
    // constexpr double alpha4 = 0.0193;

    constexpr double w_cutoff = 0.5;
    constexpr double gain = -1.0;
    constexpr double eps = 1e-6;

    // -----------------------------
    // Extract states
    // -----------------------------
    const double u = state[3];
    const double v = state[4];
    const double r = state[5];
    const double psi = state[2];
    const double F_eff = state[7];
    const double delta = state[6]*0.01;

    const double ur = u + stream_speed*cos(psi);
    const double vr = v - stream_speed*sin(psi);

    // std::cout << "delta=" << delta << "\n";
    const double T = F_eff;
    // -----------------------------
    // Nominal dynamics f(x)
    // -----------------------------
    std::array<double,3> f_usv {{
        (Xu * ur + Xuu * std::sqrt(ur*ur + eps) * ur + u_F * T * std::cos(delta))/(M + Xu_dot),

        (Yv * vr + Yr * r + Yvv * std::sqrt(vr*vr + eps) * vr + Yrr * std::sqrt(r*r + eps) * r  + Yrrr*r*r*r + v_F * T * std::sin(delta))/(M + Yv_dot),

        (Nr * r + Nv * vr + Nrr * std::sqrt(r*r + eps) * r + Nrrr*r*r*r - r_F * T * std::sin(delta))/(I + Nr_dot)
    }};

    // -----------------------------
    // Estimation error
    // -----------------------------
    std::array<double,3> x_error {{
        dob.state_estim[0] - u,
        dob.state_estim[1] - v,
        dob.state_estim[2] - r
    }};

    // -----------------------------
    // Observer dynamics
    // -----------------------------
    std::array<double,3> xdot {{
        dob.param_estim[0] + f_usv[0] + gain * x_error[0],
        dob.param_estim[1] + f_usv[1] + gain * x_error[1],
        dob.param_estim[2] + f_usv[2] + gain * x_error[2]
    }};

    // Integrate observer state
    for (int i = 0; i < 3; ++i)
        dob.state_estim[i] += xdot[i] * dt;

    // -----------------------------
    // Disturbance estimation
    // -----------------------------
    const double pi = (1.0 / gain) * (std::exp(gain * dt) - 1.0);

    for (int i = 0; i < 3; ++i)
        dob.param_estim[i] = -std::exp(gain * dt) * x_error[i] / pi;

    // -----------------------------
    // Low-pass filtering
    // -----------------------------
    const double alpha = std::exp(-w_cutoff * dt);

    for (int i = 0; i < 3; ++i)
        dob.param_filtered[i] =
            dob.param_filtered[i] * alpha
          - dob.param_estim[i] * (1.0 - alpha);
}

} // namespace DOB
