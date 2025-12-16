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
    double dt
)
{
    // -----------------------------
    // Model parameters
    // -----------------------------
    constexpr double Xu = 0.0845;
    constexpr double Xuu = 0.0195;
    constexpr double Yv = 0.0485;
    constexpr double Yvv = 0.0988;
    constexpr double Yr = 0.151;
    constexpr double Nr = 0.6939;
    constexpr double Nrr = 0.0;
    constexpr double Nv = 0.0;

    constexpr double alpha1 = 0.0452;
    constexpr double alpha2 = 2.0 / 500.0;
    constexpr double alpha3 = 0.0188;
    constexpr double alpha4 = 0.0193;

    constexpr double w_cutoff = 1.0;
    constexpr double gain = -1.0;
    constexpr double eps = 1e-6;

    // -----------------------------
    // Extract states
    // -----------------------------
    const double u = state[3];
    const double v = state[4];
    const double r = state[5];
    const double F_eff = state[7];
    const double delta = state[6];

    // -----------------------------
    // Thrust model
    // -----------------------------
    constexpr double s = 25.0;
    constexpr double k = 8.0;
    constexpr double a1 = 2.2 * 2.2;
    constexpr double a2 = 2.2 * 2.2;
    constexpr double b11 = 1.0;
    constexpr double b22 = 1.0;

    const double T =
        (1.0 / (1.0 + std::exp(s * F_eff))) *
            (b11 * F_eff + std::tanh(k * F_eff) * a1)
      + (1.0 / (1.0 + std::exp(-s * F_eff))) *
            (b22 * F_eff + std::tanh(k * F_eff) * a2);

    // -----------------------------
    // Nominal dynamics f(x)
    // -----------------------------
    std::array<double,3> f_usv {{
        -Xu * u - Xuu * std::sqrt(u*u + eps) * u
            + alpha1 * T * std::cos(alpha2 * delta),

        -Yv * v - Yr * r - Yvv * std::sqrt(v*v + eps) * v
            + alpha3 * T * std::sin(alpha2 * delta),

        -Nr * r - Nv * v - Nrr * std::sqrt(r*r + eps) * r
            - alpha4 * T * std::sin(alpha2 * delta)
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
