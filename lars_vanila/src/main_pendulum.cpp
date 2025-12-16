#include <iostream>
#include <cmath>

#include "acados/utils/print.h"
#include "acados_c/external_function_interface.h"

#include "acados_solver_pendulum.h"

int main()
{
    pendulum_solver_capsule* capsule = pendulum_acados_create_capsule();
    if (!capsule)
    {
        std::cerr << "Error: create capsule failed\n";
        return 1;
    }

    if (pendulum_acados_create(capsule))
    {
        std::cerr << "Error: create solver failed\n";
        return 1;
    }

    const int N  = PENDULUM_N;
    const int nx = PENDULUM_NX;
    const int nu = PENDULUM_NU;

    std::cout << "Created MPC: N=" << N
              << ", nx=" << nx
              << ", nu=" << nu << std::endl;

    // ---- Set initial state ----
    double x0[PENDULUM_NX] = {0.0, 0.0, 0.2, 0.0};

    ocp_nlp_in *nlp_in = capsule->nlp_in;

    ocp_nlp_constraints_model_set(
        capsule->nlp_config,
        capsule->nlp_dims,
        nlp_in,
        capsule->nlp_out,   // <-- REQUIRED
        0,
        "lbx",
        x0);
    
    ocp_nlp_constraints_model_set(
        capsule->nlp_config,
        capsule->nlp_dims,
        nlp_in,
        capsule->nlp_out,   // <-- REQUIRED
        0,
        "ubx",
        x0);

    // ---- Solve ----
    int status = pendulum_acados_solve(capsule);
    std::cout << "Solve status: " << status << std::endl;

    if (status)
    {
        pendulum_acados_free(capsule);
        pendulum_acados_free_capsule(capsule);
        return 1;
    }

    ocp_nlp_out *nlp_out = capsule->nlp_out;

    // ---- get first control and state ----
    double u0[PENDULUM_NU];
    double x0_sol[PENDULUM_NX];

    ocp_nlp_out_get(
        capsule->nlp_config,
        capsule->nlp_dims,
        nlp_out,
        0,
        "u",
        u0);

    ocp_nlp_out_get(
        capsule->nlp_config,
        capsule->nlp_dims,
        nlp_out,
        0,
        "x",
        x0_sol);

    std::cout << "u(0): ";
    for (int i=0;i<nu;i++) std::cout << u0[i] << " ";
    std::cout << "\n";

    std::cout << "x(0): ";
    for (int i=0;i<nx;i++) std::cout << x0_sol[i] << " ";
    std::cout << "\n";

    // -------------------------------------
    // Print predicted trajectory (x(k), u(k))
    // -------------------------------------
    std::cout << "\n=== Predicted MPC Trajectory ===\n";

    for (int stage = 0; stage <= N; stage++)
    {
        double xk[PENDULUM_NX];
        ocp_nlp_out_get(
            capsule->nlp_config,
            capsule->nlp_dims,
            nlp_out,
            stage,
            "x",
            xk);

        std::cout << "x[" << stage << "] = [ ";
        for (int i = 0; i < nx; i++)
            std::cout << xk[i] << " ";
        std::cout << "]";

        if (stage < N)
        {
            double uk[PENDULUM_NU];
            ocp_nlp_out_get(
                capsule->nlp_config,
                capsule->nlp_dims,
                nlp_out,
                stage,
                "u",
                uk);

            std::cout << ",   u[" << stage << "] = [ ";
            for (int i = 0; i < nu; i++)
                std::cout << uk[i] << " ";
            std::cout << "]";
        }

        std::cout << "\n";
    }

    std::cout << "===================================\n";

    // ---- Cleanup ----
    pendulum_acados_free(capsule);
    pendulum_acados_free_capsule(capsule);

    return 0;
}
