// mpc.hpp
// - Pure MPC module (NO ROS2)
// - Owns ACADOS solver capsule
// - Loads parameters from a key=value txt (same file can be used by DOB)
// - Exposes setters for measured/estimated state and (optional) LARS info
// - step() produces updated (delta, F)

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "reference_generator.hpp"
#include "acados_solver_heron.h"

// ================================
// NOTE
// - This header intentionally contains NO ROS2 dependencies.
// - Any logging should be handled by your ROS2 main node.
// ================================

class MPCModule
{
public:

    double lars_mode() const { return lars_mode_; }

    struct Weights
    {
        // Stage Q
        double Q_x=1e3, Q_y=1e3, Q_psi=1e5, Q_u=1e1, Q_v=1e3, Q_r=1e2, Q_delta=1e1, Q_F=1e-1;
        // Terminal QN
        double QN_x=1e5, QN_y=1e5, QN_psi=1e4, QN_u=1e-2, QN_v=1e-1, QN_r=1e2, QN_delta=1e1, QN_F=1e-1;
        // Input R
        double R_ddelta=1e2, R_dF=1e-2;
    };

    struct Model
    {
        double Xu=0.0, Xuu=-1106.0;
        double Yv=-296.83, Yvv=0.0;
        double Yr=-123.0, Yrr=0.0, Yrrr=-231.0;
        double Nr=-3611.0, Nrr=0.0, Nrrr=0.0;
        double Nv=-3809.0;
        double u_F=140.0, v_F=140.0, r_F=490.0;
        double Xu_dot=1041.0, Yv_dot=749.0, Nr_dot=1066.0;
        double M=4282.0, I=9050.0;
    };

    struct Misc
    {
        // offsets
        double offset_x=290584.355;
        double offset_y=4117794.41;

        // docking point rel (in CD frame)
        double DP_x_rel=-12.0;
        double DP_y_rel=-1.0;
        double DP_psi_rel=0.0;

        // cradle / docking
        double cradle_x_rel=4.0;
        double dock_radius=2.0;
        double duration_time=5.0;
        double dp_dist=1.5;

        // current / mother ship speed
        double stream_speed=0.0;
        double MS_speed=0.0;

        // end phase
        double final_boost_radius=3.0;
        double final_end_radius=1.0;

        // hook point
        double hp=3.0;

        // CBF / misc params
        double alpha_side=0.1;
        double alpha_dock=0.5;


        // a = 1.0     # hopper half-height
        // b = 0.2     # 곡률 (얼마나 급하게 열리는지)
        // x0 = 0.0    # hopper 목 위치
        // y0 = 0.7  
        double cbf_d=0.0;
        double cbf_a=1.0;
        double cbf_b=0.2;
        double cbf_x0=-0.0;
        double cbf_y0=0.7;

        double wpt_x1 = 0.0;
        double wpt_x2 = 0.0;
        double wpt_x3 = 0.0;
        double wpt_x4 = 0.0;
        
        double wpt_y1 = 0.0;
        double wpt_y2 = 0.0;
        double wpt_y3 = 0.0;
        double wpt_y4 = 0.0;
    
        double acceptance_radius = 0.5;
    };

    struct Params
    {
        Weights w;
        Model   m;
        Misc    misc;
    };

    struct LarsState
    {
        double CD_x=0.0, CD_y=0.0, CD_psi=0.0;
        double MS_x=0.0, MS_y=0.0, MS_psi=0.0;
    };

    // ================================
    // Ctor / Dtor
    // ================================
    MPCModule() = default;

    ~MPCModule()
    {
        if (capsule_)
        {
            heron_acados_free(capsule_);
            heron_acados_free_capsule(capsule_);
            capsule_ = nullptr;
        }
    }

    // Call once after construction.
    // - param_txt_path: key=value file (weights/model/misc)
    // - brs_csv_path: reachability grid file
    // Returns true on success.
    bool init(
        const std::string& param_txt_path,
        const std::string& brs_csv_path,
        int N = HERON_N,
        double dt_control = 0.25)
    {
        param_txt_path_ = param_txt_path;
        N_ = N;
        dt_control_ = dt_control;

        // Create ACADOS capsule
        capsule_ = heron_acados_create_capsule();
        if (!capsule_) return false;
        if (heron_acados_create(capsule_)) return false;

        // BRS grid
        brs_.dx   = (brs_.xmax   - brs_.xmin)   / (brs_.Nx - 1);
        brs_.dy   = (brs_.ymax   - brs_.ymin)   / (brs_.Nx - 1);
        brs_.dpsi = (brs_.psimax - brs_.psimin) / (brs_.Nx - 1);

        if (!loadBRSFromCSV(brs_csv_path)) return false;

        // Load params and apply cost
        (void)loadParamsFromTxt(param_txt_path_);
        applyCostWeights();

        // Reference generator defaults (kept from your original)
        A_ = 6;
        B_ = 3;
        C_ = 6;
        theta_ = 0.0;
        ref_dt_ = 0.001;

        // Initialize offsets if not provided
        // (offsets are used as origin for relative states)
        // They can be overwritten by txt.

        // Warm-up flag
        warmed_up_ = false;

        return true;
    }

    // ================================
    // Public setters
    // ================================
    void setEstimatedState(double x, double y, double psi, double u, double v, double r)
    {
        x_ = x;
        y_ = y;
        yaw_ = psi;
        u_ = u;
        v_ = v;
        r_ = r;


    }

    // mpc.hpp (inside class MPCModule, public:)
    bool getPredictedXY(std::vector<double>& out_x,
        std::vector<double>& out_y,
        bool absolute = true) const
    {
    if (!capsule_) return false;

    out_x.resize((size_t)N_ + 1);
    out_y.resize((size_t)N_ + 1);

    for (int j = 0; j <= N_; ++j)
    {
    double xj[HERON_NX];
    ocp_nlp_out_get(
    capsule_->nlp_config,
    capsule_->nlp_dims,
    capsule_->nlp_out,
    j,
    "x",
    xj);

    // xj[0], xj[1]은 (상대좌표) 상태이므로 absolute면 offset을 더해 절대좌표로
    out_x[(size_t)j] = xj[0] + (absolute ? params_.misc.offset_x : 0.0);
    out_y[(size_t)j] = xj[1] + (absolute ? params_.misc.offset_y : 0.0);
    }


    ////////////// reference trajectory //////////////
    // for (int k = 0; k <= N_; ++k)
    // {

    //     int j = ref_idx_ + k * ref_step_;
    //     // ref_idx_ += 1;
    //     double xr = reference_[j].x;
    //     double yr = reference_[j].y;

    //     // reference_가 "상대좌표"라면 absolute일 때 offset 더하기
    

    //     out_x[static_cast<size_t>(k)] = xr;
    //     out_y[static_cast<size_t>(k)] = yr;
    //     // std::cout << "xr :" << xr - x_ + params_.misc.offset_x << "yr :" << yr - y_ + params_.misc.offset_y  <<  std::endl;
    //     // std::cout << "xr :" << xr - params_.misc.offset_x<< ", yr :" << yr - params_.misc.offset_y  <<  std::endl;
    // }
    // std::cout << "xr :" << out_x[0]<< " 5 :" << out_x[15]<< " 10 :" << out_x[30] <<std::endl;


    return true;
    }

    void setActuatorInternal(double delta, double F)
    {
        delta_ = delta;
        F_ = F;
    }

    void setLarsState(const LarsState& ls)
    {
        lars_ = ls;
    }



    double dt_control() const { return dt_control_; }
    int N() const { return N_; }

    // Expose current integrated actuator states (needed by DOB + ROS publishing)
    double delta() const { return delta_; }
    double F() const { return F_; }

    const Params& params() const { return params_; }
    Params& params() { return params_; }

    // ================================
    // Main step
    // ================================
    // Inputs:
    // - dob_u/dob_v/dob_r: disturbance estimates (filtered)
    // - now_sec: current time in seconds (monotonic or ROS time ok)
    // Outputs:
    // - delta_out, F_out : updated integrated actuator states
    // Return true if solver ran.
    bool step(double dob_u, double dob_v, double dob_r, double now_sec, double& delta_out, double& F_out)
    {
        // std::fprintf(stderr,
        //     "[Lars mode -1] =%.3f\n",
        //     lars_mode_);

        if (!capsule_) return false;

        // Periodic reload of txt (every 5 seconds, like your original)
        if (!time_init_)
        {
            time_init_ = true;
            start_time_sec_ = now_sec;
            last_param_reload_sec_ = now_sec;
        }
        if ((now_sec - last_param_reload_sec_) >= 5.0)
        {
            (void)loadParamsFromTxt(param_txt_path_);
            applyCostWeights();
            last_param_reload_sec_ = now_sec;
        }

        // Build relative state vector for ACADOS (x,y are relative to offset)
        updateStateVector();

        // 1) WARM-UP: do once
        if (!warmed_up_)
        {
            // set ref = current
            ref_x_ = x_;
            ref_y_ = y_;
            theta_ = yaw_;

            setInitialStateBounds();

            int rti_phase = 1;
            ocp_nlp_solver_opts_set(capsule_->nlp_config, capsule_->nlp_opts, "rti_phase", &rti_phase);
            heron_acados_solve(capsule_);

            generateReference();
            warmed_up_ = true;
            delta_out = delta_;
            F_out = F_;
            return true;
        }

        // 2) Reachability + mode check (kept; relies on LARS inputs)



        BRS_exist_ = updateBRSExist();


        ModeCheck();


        checkWaypointHoldAndSwitch();


        // 3) Reference (kept from original)
        // - If eight_traj_check==false: straight waypoint following
        // - Else: figure-eight
        for (int j = 0; j <= N_; ++j)
        {
            double xr = 0.0;
            double yr = 0.0;
            double pr = 0.0;

            if (!eight_traj_check_)
            {
                xr = wpt_x_;
                yr = wpt_y_;
                pr = wpt_psi_;
            }
            else
            {
                if (!reference_.empty())
                {
                    int k = std::min(ref_idx_ + j * ref_step_, (int)reference_.size() - 1);

                
                    xr = params_.misc.DP_x_rel;
                    yr = params_.misc.DP_y_rel;
                    pr = 0.0;
                }
            }

            if (j == N_)
            {
                double yrefN[HERON_NYN] = {xr, yr, pr, 0.0, 0,0,0,0};
                ocp_nlp_cost_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, j, "yref", yrefN);
            }
            else
            {
                double yref[HERON_NY] = {xr, yr, pr, 0.0, 0,0,0,0,0,0};
                ocp_nlp_cost_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, j, "yref", yref);
            }
        }
        ref_idx_ = ref_idx_ + 1;

        // 4) Prepare parameter vector (same ordering as original)
        // Apply DOB after 5 seconds from start
        const double t_from_start = now_sec - start_time_sec_;
        const double du = (t_from_start >= 5.0) ? dob_u : 0.0;
        const double dv = (t_from_start >= 5.0) ? dob_v : 0.0;
        const double dr = (t_from_start >= 5.0) ? dob_r : 0.0;

        const double stream_speed = params_.misc.stream_speed;
        const double switch_1 = (lars_mode_ >= 1.0) ? 1.0 : 0.0;

        lars_.CD_x = params_.misc.wpt_x1;
        lars_.CD_y = params_.misc.wpt_y1;
        lars_.CD_psi = 0.0;
        
        lars_.MS_x = params_.misc.wpt_x2;
        lars_.MS_y = params_.misc.wpt_y2;
        lars_.MS_psi = 0.0;
        

        double obs[HERON_NP] = {
            lars_.CD_x - params_.misc.offset_x, lars_.CD_y - params_.misc.offset_y, lars_.CD_psi,
            lars_.MS_x - params_.misc.offset_x, lars_.MS_y - params_.misc.offset_y, lars_.MS_psi,
            du, dv, dr,
            // 0.0,0.0,0.0,
            stream_speed, switch_1,
            params_.m.Xu, params_.m.Xuu,
            params_.m.Yv, params_.m.Yvv, params_.m.Yr, params_.m.Yrr, params_.m.Yrrr,
            params_.m.Nr, params_.m.Nrr, params_.m.Nrrr, params_.m.Nv,
            params_.m.u_F, params_.m.v_F, params_.m.r_F,
            params_.m.Xu_dot, params_.m.Yv_dot, params_.m.Nr_dot,
            params_.m.M, params_.m.I,
            params_.misc.alpha_side, params_.misc.alpha_dock,
            params_.misc.cbf_d, params_.misc.cbf_a, params_.misc.cbf_b, params_.misc.cbf_x0, params_.misc.cbf_y0
        };



        for (int j = 0; j <= N_; j++)
            heron_acados_update_params(capsule_, j, obs, HERON_NP);

        // 5) State equality constraint update
        setInitialStateBounds();

        // 6) RTI solve
        int rti_phase = 1;
        ocp_nlp_solver_opts_set(capsule_->nlp_config, capsule_->nlp_opts, "rti_phase", &rti_phase);
        heron_acados_solve(capsule_);

        rti_phase = 2;
        ocp_nlp_solver_opts_set(capsule_->nlp_config, capsule_->nlp_opts, "rti_phase", &rti_phase);
        heron_acados_solve(capsule_);

        // 7) Get control output
        double uMPC[HERON_NU];
        ocp_nlp_out_get(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_out, 0, "u", uMPC);

        delta_ += uMPC[0] * dt_control_;
        F_     += uMPC[1] * dt_control_;

        // clamp (kept)
        if (F_ >= 3.5) F_ = 3.5 - 0.1;
        else if (F_ <= -1.3) F_ = -1.5 + 0.01;

        const double delta_max = 3000.0 * M_PI / 180.0;
        if (delta_ >  delta_max) delta_ =  delta_max - 0.2;
        if (delta_ < -delta_max) delta_ = -delta_max + 0.2;


        std::cout << "WPT x : "<< wpt_x_ << ", y : " << wpt_y_ <<", Distance : "<<std::sqrt((wpt_x_ - x_)*(wpt_x_ - x_) + (wpt_y_ - y_)*(wpt_y_ - y_))<<std::endl;
        std::cout << "Lars mode : "<< lars_mode_ << std::endl;

        // if std::sqrt()

        delta_out = delta_;
        F_out = F_;
        return true;
    }

private:
    // ================================
    // TXT key=value parser (local)
    // ================================

    double lars_mode_ = 0.0;   // mode는 여기서만 관리

    static inline std::string trim_copy(std::string s)
    {
        auto is_space = [](unsigned char c){ return std::isspace(c) != 0; };
        while (!s.empty() && is_space(s.front())) s.erase(s.begin());
        while (!s.empty() && is_space(s.back()))  s.pop_back();
        return s;
    }

    bool loadParamsFromTxt(const std::string& path)
    {
        std::ifstream file(path);
        if (!file.is_open()) return false;

        std::string line;
        while (std::getline(file, line))
        {
            line = trim_copy(line);
            if (line.empty() || line[0] == '#') continue;

            auto pos = line.find('=');
            if (pos == std::string::npos) continue;

            std::string key = trim_copy(line.substr(0, pos));
            std::string val = trim_copy(line.substr(pos + 1));

            try
            {
                double v = std::stod(val);

                // misc
                if (key == "DP_x_rel") params_.misc.DP_x_rel = v;
                else if (key == "DP_y_rel") params_.misc.DP_y_rel = v;
                else if (key == "DP_psi_rel") params_.misc.DP_psi_rel = v;
                else if (key == "offset_x") params_.misc.offset_x = v;
                else if (key == "offset_y") params_.misc.offset_y = v;
                else if (key == "cradle_x_rel") params_.misc.cradle_x_rel = v;
                else if (key == "duration_time") params_.misc.duration_time = v;
                else if (key == "dp_dist") params_.misc.dp_dist = v;
                else if (key == "stream_speed") params_.misc.stream_speed = v;
                else if (key == "dock_radius") params_.misc.dock_radius = v;
                else if (key == "MS_speed") params_.misc.MS_speed = v;
                else if (key == "alpha_side") params_.misc.alpha_side = v;
                else if (key == "alpha_dock") params_.misc.alpha_dock = v;
                else if (key == "cbf_d") params_.misc.cbf_d = v;
                else if (key == "cbf_a") params_.misc.cbf_a = v;
                else if (key == "cbf_b") params_.misc.cbf_b = v;
                else if (key == "cbf_x0") params_.misc.cbf_x0 = v;
                else if (key == "cbf_y0") params_.misc.cbf_y0 = v;
                else if (key == "hp") params_.misc.hp = v;
                else if (key == "final_boost_radius") params_.misc.final_boost_radius = v;
                else if (key == "final_end_radius") params_.misc.final_end_radius = v;

                // model
                else if (key == "Xu") params_.m.Xu = v;
                else if (key == "Xuu") params_.m.Xuu = v;
                else if (key == "Yv") params_.m.Yv = v;
                else if (key == "Yvv") params_.m.Yvv = v;
                else if (key == "Yr") params_.m.Yr = v;
                else if (key == "Yrr") params_.m.Yrr = v;
                else if (key == "Yrrr") params_.m.Yrrr = v;
                else if (key == "Nr") params_.m.Nr = v;
                else if (key == "Nrr") params_.m.Nrr = v;
                else if (key == "Nrrr") params_.m.Nrrr = v;
                else if (key == "Nv") params_.m.Nv = v;
                else if (key == "u_F") params_.m.u_F = v;
                else if (key == "v_F") params_.m.v_F = v;
                else if (key == "r_F") params_.m.r_F = v;
                else if (key == "Xu_dot") params_.m.Xu_dot = v;
                else if (key == "Yv_dot") params_.m.Yv_dot = v;
                else if (key == "Nr_dot") params_.m.Nr_dot = v;
                else if (key == "M") params_.m.M = v;
                else if (key == "I") params_.m.I = v;

                // weights
                else if (key == "Q_x") params_.w.Q_x = v;
                else if (key == "Q_y") params_.w.Q_y = v;
                else if (key == "Q_psi") params_.w.Q_psi = v;
                else if (key == "Q_u") params_.w.Q_u = v;
                else if (key == "Q_v") params_.w.Q_v = v;
                else if (key == "Q_r") params_.w.Q_r = v;
                else if (key == "Q_delta") params_.w.Q_delta = v;
                else if (key == "Q_F") params_.w.Q_F = v;

                else if (key == "QN_x") params_.w.QN_x = v;
                else if (key == "QN_y") params_.w.QN_y = v;
                else if (key == "QN_psi") params_.w.QN_psi = v;
                else if (key == "QN_u") params_.w.QN_u = v;
                else if (key == "QN_v") params_.w.QN_v = v;
                else if (key == "QN_r") params_.w.QN_r = v;
                else if (key == "QN_delta") params_.w.QN_delta = v;
                else if (key == "QN_F") params_.w.QN_F = v;

                else if (key == "R_ddelta") params_.w.R_ddelta = v;
                else if (key == "R_dF") params_.w.R_dF = v;

                else if (key == "wpt_x1") params_.misc.wpt_x1 = v;
                else if (key == "wpt_x2") params_.misc.wpt_x2 = v;
                else if (key == "wpt_x3") params_.misc.wpt_x3 = v;
                else if (key == "wpt_x4") params_.misc.wpt_x4 = v;


                else if (key == "wpt_y1") params_.misc.wpt_y1 = v;
                else if (key == "wpt_y2") params_.misc.wpt_y2 = v;
                else if (key == "wpt_y3") params_.misc.wpt_y3 = v;
                else if (key == "wpt_y4") params_.misc.wpt_y4 = v;

                else if (key == "acceptance_radius") params_.misc.acceptance_radius;


                // ignore others
            }
            catch (...) { /* ignore */ }
        }

        return true;
    }

    void applyCostWeights()
    {
        // Stage W (HERON_NY x HERON_NY)
        std::vector<double> W(HERON_NY * HERON_NY, 0.0);
        auto set_diag = [&](int idx, double val){
            W[idx * HERON_NY + idx] = 2.0 * val;
        };

        set_diag(0, params_.w.Q_x);
        set_diag(1, params_.w.Q_y);
        set_diag(2, params_.w.Q_psi);
        set_diag(3, params_.w.Q_u);
        set_diag(4, params_.w.Q_v);
        set_diag(5, params_.w.Q_r);
        set_diag(6, params_.w.Q_delta);
        set_diag(7, params_.w.Q_F);
        set_diag(8, params_.w.R_ddelta);
        set_diag(9, params_.w.R_dF);

        for (int j = 0; j < N_; ++j)
        {
            ocp_nlp_cost_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, j, "W", W.data());
        }

        // Terminal WN (HERON_NYN x HERON_NYN)
        std::vector<double> WN(HERON_NYN * HERON_NYN, 0.0);
        auto set_diagN = [&](int idx, double val){
            WN[idx * HERON_NYN + idx] = 2.0 * val;
        };

        set_diagN(0, params_.w.QN_x);
        set_diagN(1, params_.w.QN_y);
        set_diagN(2, params_.w.QN_psi);
        set_diagN(3, params_.w.QN_u);
        set_diagN(4, params_.w.QN_v);
        set_diagN(5, params_.w.QN_r);
        set_diagN(6, params_.w.QN_delta);
        set_diagN(7, params_.w.QN_F);

        ocp_nlp_cost_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, N_, "W", WN.data());
    }

    // =====================================
    // BRS
    // =====================================
    struct BRSGrid {
        int Nx = 21;
        std::vector<double> data;   // size = Nx^4
        double xmin = -10.0, xmax = 2.0;
        double ymin = -6.0, ymax = 6.0;
        double psimin = -3.141592, psimax = 3.141592;
        double dx = 0.0, dy = 0.0, dpsi = 0.0;
    };

    bool loadBRSFromCSV(const std::string& path)
    {
        std::ifstream file(path);
        if (!file.is_open()) return false;

        brs_.data.clear();
        brs_.data.reserve((size_t)std::pow(brs_.Nx, 4));

        std::string line;
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ','))
            {
                try { brs_.data.push_back(std::stod(cell)); }
                catch (...) {}
            }
        }
        file.close();

        const size_t expected = (size_t)std::pow(brs_.Nx, 4);
        return brs_.data.size() == expected;
    }

    bool updateBRSExist()
    {
        // 1) hook point in world
        const double px = x_ + params_.misc.hp * std::cos(yaw_);
        const double py = y_ + params_.misc.hp * std::sin(yaw_);

        // 2) cradle point in world
        const double cx = (lars_.CD_x) + params_.misc.cradle_x_rel * std::cos(lars_.CD_psi);
        const double cy = (lars_.CD_y) + params_.misc.cradle_x_rel * std::sin(lars_.CD_psi);

        // 3) relative vector in world
        const double dx = px - cx;
        const double dy = py - cy;

        // 4) world -> CD frame
        const double c = std::cos(lars_.CD_psi);
        const double s = std::sin(lars_.CD_psi);
        const double x_rel =  c * dx + s * dy;
        const double y_rel = -s * dx + c * dy;

        if (std::sqrt(x_rel*x_rel + y_rel*y_rel) < params_.misc.dock_radius || x_rel > params_.misc.dock_radius)
            return true;

        // 2) bounds
        if (x_rel < brs_.xmin || x_rel > brs_.xmax || y_rel < brs_.ymin || y_rel > brs_.ymax)
            return false;

        // 3) index
        auto idx = [&](double x, double xmin, double dx) {
            return static_cast<int>(std::round((x - xmin) / dx));
        };
        const int ix = idx(x_rel, brs_.xmin, brs_.dx);
        const int iy = idx(y_rel, brs_.ymin, brs_.dy);
        if (ix < 0 || ix >= brs_.Nx || iy < 0 || iy >= brs_.Nx) return false;

        // 4) lookup min over psi,w dims
        auto brsIndex = [&](int ix, int iy, int iz, int iw) {
            return (((iw * brs_.Nx + iz) * brs_.Nx + iy) * brs_.Nx + ix);
        };

        double phi_min = std::numeric_limits<double>::infinity();
        for (int iz = 0; iz < brs_.Nx; ++iz)
            for (int iw = 0; iw < brs_.Nx; ++iw)
                phi_min = std::min(phi_min, brs_.data[(size_t)brsIndex(ix, iy, iz, iw)]);

        // std::fprintf(stderr, "Inside BRS = ", (phi_min < 0.0) );
        // std::fprintf(stderr, "Inside BRS = %s\n", (phi_min < 0.0) ? "TRUE" : "FALSE");

        return (phi_min < 0.0);
    }

    // =====================================
    // Reference
    // =====================================
    void generateReference()
    {
        reference_ = ReferenceGenerator::generateFigureEight(10000.0, ref_dt_, ref_x_, ref_y_, theta_, A_, B_, C_);
        std::vector<RefPoint> reduced;
        int step = (int)(dt_control_ / ref_dt_);
        step = std::max(1, step);
        for (int i = 0; i < (int)reference_.size(); i += step)
            reduced.push_back(reference_[i]);
        reference_ = reduced;
    }

    // =====================================
    // State vectors for constraints
    // =====================================
    void updateStateVector()
    {
        states_[0] = x_ - params_.misc.offset_x;
        states_[1] = y_ - params_.misc.offset_y;
        states_[2] = yaw_;
        states_[3] = u_;
        states_[4] = v_;
        states_[5] = r_;
        states_[6] = delta_;
        states_[7] = F_;

        // std::fprintf(stderr,
        //     "[MPC setEstimatedState] x=%.3f y=%.3f psi=%.3f u=%.3f v=%.3f r=%.3f\n",
        //     states_[0], states_[1], states_[2], states_[3], states_[4], states_[5]);
    }

    void setInitialStateBounds()
    {
        ocp_nlp_constraints_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, capsule_->nlp_out, 0, "lbx", states_);
        ocp_nlp_constraints_model_set(capsule_->nlp_config, capsule_->nlp_dims, capsule_->nlp_in, capsule_->nlp_out, 0, "ubx", states_);
    }



    // =====================================
    // Mode / waypoint logic (kept)
    // =====================================
    void ModeCheck()
    {

        // std::fprintf(stderr,
        //     "[Lars mode] =%.3f\n",
        //     lars_.lars_mode);

        // BRS out -> return to mode 0
        // if (!BRS_exist_ && lars_mode_ == 1.0)
        // {
        //     // lars_mode_ = 0.0;
        //     lars_auto_triggered_ = false;
        //     inside_wpt_radius_ = false;
            
        //     wpt_x_ = params_.misc.DP_x_rel;
        //     wpt_y_ = params_.misc.DP_y_rel;
        //     wpt_psi_ = 0.0;

        //     return;
        // }

        if (lars_mode_ == 0.0)
        {
            wpt_x_ = params_.misc.DP_x_rel ;
            wpt_y_ = params_.misc.DP_y_rel ;
            wpt_psi_ = 0.0;

        }
        else if (lars_mode_ == 1.0)
        {
            wpt_x_ = params_.misc.wpt_x1;
            wpt_y_ = params_.misc.wpt_y1;
            wpt_psi_ = 0.0;

        }
        else if (lars_mode_ == 2.0)
        {
            wpt_x_ = params_.misc.wpt_x1;
            wpt_y_ = params_.misc.wpt_y1;
            wpt_psi_ = 0.0;

        }
    }

    void checkWaypointHoldAndSwitch()
    {
        // terminal acting
        const double hp_x = x_ + params_.misc.hp * std::cos(yaw_) - params_.misc.offset_x;
        const double hp_y = y_ + params_.misc.hp * std::sin(yaw_) - params_.misc.offset_y;

        // if (std::sqrt((hp_x - lars_.CD_x + params_.misc.offset_x)*(hp_x - lars_.CD_x + params_.misc.offset_x) + (hp_y - lars_.CD_y + params_.misc.offset_y)*(hp_y - lars_.CD_y + params_.misc.offset_y)) < params_.misc.final_boost_radius)
        // {
        //     lars_mode_ = 2.0;
        //     final_get_in_ = true;
        // }

        // // if (lars_auto_triggered_) return;
        // if (!BRS_exist_) { inside_wpt_radius_ = false; return; }

        // distance to waypoint in global frame
        const double dx = x_ - (wpt_x_ + params_.misc.offset_x);
        const double dy = y_ - (wpt_y_ + params_.misc.offset_y);
        const double dist = std::sqrt(dx*dx + dy*dy);



        // NOTE: Since this module has no ROS time, we approximate the hold timer
        // using the control dt. This is good enough for switching logic.
        if ((dist <= params_.misc.dp_dist) && (lars_mode_ == 0.0))
        {
            if (!inside_wpt_radius_)
            {
                inside_wpt_radius_ = true;
                wpt_hold_elapsed_ = 0.0;
            }
            else
            {
                wpt_hold_elapsed_ += dt_control_;
                if (wpt_hold_elapsed_ >= params_.misc.duration_time)
                {
                    lars_mode_ = 1.0;
                    lars_auto_triggered_ = true;
                    std::fprintf(stderr, "[Start docking] \n");
                }
            }
        }
        else
        {
            inside_wpt_radius_ = false;
            wpt_hold_elapsed_ = 0.0;
        }

        std::fprintf(stderr, "[dist] = %f, [time] = %f \n", dist, wpt_hold_elapsed_);
    }

private:
    // solver
    heron_solver_capsule* capsule_ = nullptr;

    // configuration
    std::string param_txt_path_;
    int N_ = HERON_N;
    double dt_control_ = 0.25;
    Params params_{};
    LarsState lars_{};

    // measured/estimated state
    double x_ = 0.0, y_ = 0.0, yaw_ = 0.0;
    double u_ = 0.0, v_ = 0.0, r_ = 0.0;
    double delta_ = 0.0, F_ = 0.0;

    // constraint vector
    double states_[HERON_NX] = {0};

    // reference
    std::vector<RefPoint> reference_;
    double A_=80, B_=80, C_=1000, theta_=0.0;
    double ref_dt_=0.001;
    double ref_x_=0.0, ref_y_=0.0;

    bool eight_traj_check_ = false;
    int ref_idx_ = 0;
    int ref_step_ = 1;

    // warm-up
    bool warmed_up_ = false;

    // timing
    bool time_init_ = false;
    double start_time_sec_ = 0.0;
    double last_param_reload_sec_ = 0.0;

    // BRS
    BRSGrid brs_{};
    bool BRS_exist_ = false;

    // switching
    bool inside_wpt_radius_ = false;
    bool lars_auto_triggered_ = false;
    double wpt_hold_elapsed_ = 0.0;

    bool final_get_in_ = false;
    bool mission_over_ = false;

    // waypoint
    double wpt_x_ = 0.0;
    double wpt_y_ = -13.0;
    double wpt_psi_ = 0.0;


};
