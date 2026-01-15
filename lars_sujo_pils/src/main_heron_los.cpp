#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <std_msgs/msg/float64.hpp>

#include "reference_generator.hpp"
#include "dob.hpp"
#include "uvr_filter.hpp"
#include "acados_solver_heron.h"

#include "aura_msg/msg/mpc_state.hpp"
#include "aura_msg/msg/mpc_traj.hpp"
#include "aura_msg/msg/obs_state.hpp"

#include <cmath>
#include <chrono>
#include <vector>

#include <fstream>
#include <sstream>
#include <string>
#include <limits>


using namespace std::chrono_literals;

static constexpr double ship_state_x = (289577.66 + 291591.05) * 0.5;
static constexpr double ship_state_y = (4117065.30 + 4118523.52) * 0.5;

class AuraMPC : public rclcpp::Node
{
public:
    AuraMPC() : Node("aura_mpc")
    {
        actuator_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>(
            "/actuator_outputs", 10);

        ekf_sub_ = this->create_subscription<std_msgs::msg::Float64MultiArray>(
            "/ekf/estimated_state", 10,
            std::bind(&AuraMPC::ekf_callback, this, std::placeholders::_1));

        lars_state_sub_ = this->create_subscription<std_msgs::msg::Float64MultiArray>(
            "/lars_state", 10,
            std::bind(&AuraMPC::lars_state_callback, this, std::placeholders::_1));            

        lars_mode_sub_ = this->create_subscription<std_msgs::msg::Float64>(
            "/lars_mode", 10,
            std::bind(&AuraMPC::lars_mode_callback, this, std::placeholders::_1));  

        mpcvis_pub_ = this->create_publisher<aura_msg::msg::MPCTraj>("/mpc_vis", 10);
        dob_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>(
            "/DOB", 10);

        // Initial states
        x_ = ship_state_x;
        y_ = ship_state_y;
        yaw_ = 0.0;
        u_ = v_ = r_ = 0.0;
        delta_ = 0.0;
        F_ = 0.0;
        ref_x_ = 0.0;
        ref_y_ = 0.0;

        offset_x_ = ship_state_x;
        offset_y_ = ship_state_y;

        CD_x = ship_state_x;
        CD_y = ship_state_y;
        CD_psi = 0.0;

        MS_x = ship_state_x;
        MS_y = ship_state_y;
        MS_psi = 0.0;

        lars_mode_ = 0.0;
        switch_1 = 0.0;
        stream_speed = 0.0;

        wpt_x = 0.0;
        wpt_y = -13.0;
        wpt_psi = 0.0;

        DP_x_rel = -12.0;
        DP_y_rel = -1.0;
        DP_psi_rel = 0.0;
        
        cradle_x_rel = 4.0;
        dock_radius = 2.0;
        duration_time = 5.0;
        dp_dist = 1.5;

        N_ = HERON_N;
        dt_control_ = 0.25;

        Xu      = -0.0;
        Xuu     = -1106;
        Yv      = -296.83;
        Yvv     = -0.0;
        Yr      = -123;
        Yrr     = -0.0;
        Yrrr    = -231;
        Nr      = -3611;
        Nrr     = -0.0;
        Nrrr    = -0.0;
        Nv      = -3809;
        u_F     = 140;
        v_F     = 140;
        r_F     = 490;

        Xu_dot = 1041;
        Yv_dot = 749;
        Nr_dot = 1066;
        M = 4282;
        I = 9050;

        alpha_side = 0.1;
        alpha_dock = 0.5;
    
        cbf_d = 4;
        cbf_a = 3.0;
        cbf_b = 0.5;
        cbf_x0 = -10.0;
        cbf_y0 = 3.5;

        only_xypsi = true;

        double final_boost_radius = 3.0;
        double final_end_radius = 1.0;
    
        double hp = 3;

        MS_speed = 0.0;
    
        // === Create ACADOS solver ===
        capsule_ = heron_acados_create_capsule();

        int status = heron_acados_create(capsule_);
        if (status)
            RCLCPP_ERROR(this->get_logger(), "ACADOS solver creation FAILED!");
        else
            RCLCPP_INFO(this->get_logger(), "ACADOS solver created.");

        // reference trajectory parameters
        A_ = 120;
        B_ = 85;
        C_ = 28;
        theta_ = 0.0;
        ref_dt_ = 0.01;

        // ✅ (1) 반드시 한 번 상태벡터 초기화
        updateStateVector();

        timer_ = this->create_wall_timer(
            std::chrono::duration<double>(dt_control_),
            std::bind(&AuraMPC::run, this));

        dob_timer_ = this->create_wall_timer(
            std::chrono::duration<double>(0.1),   // 10 Hz
            std::bind(&AuraMPC::runDOB, this)
        );            

        start_time_ = this->get_clock()->now();
        last_param_update_time_ = this->get_clock()->now();

        RCLCPP_INFO(this->get_logger(), "AuraMPC initialized.");

        // ===============================
        // BRS grid initialization
        // ===============================
        brs_.dx   = (brs_.xmax   - brs_.xmin)   / (brs_.Nx - 1);
        brs_.dy   = (brs_.ymax   - brs_.ymin)   / (brs_.Nx - 1);
        brs_.dpsi = (brs_.psimax - brs_.psimin) / (brs_.Nx - 1);

        // ---- CSV path ----
        std::string brs_csv =
            "/home/user/aura_ws/src/lars/BRS_file/brs.csv";

        if (!loadBRSFromCSV(brs_csv))
        {
            RCLCPP_FATAL(this->get_logger(), "[BRS] Load failed");
            rclcpp::shutdown();
            return;
        }

        loadParamsFromTxt(param_txt_path_);
        applyCostWeights(); 

    }

    ~AuraMPC()
    {
        heron_acados_free(capsule_);
        heron_acados_free_capsule(capsule_);
    }

private:
    // ================================
    // TXT key=value parser
    // ================================
    static inline std::string trim_copy(std::string s)
    {
        auto is_space = [](unsigned char c){ return std::isspace(c); };
        while (!s.empty() && is_space(s.front())) s.erase(s.begin());
        while (!s.empty() && is_space(s.back()))  s.pop_back();
        return s;
    }

    bool loadParamsFromTxt(const std::string& path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            RCLCPP_ERROR(this->get_logger(),
                "[PARAM] Cannot open txt: %s", path.c_str());
            return false;
        }

        // 기본값: "현재 멤버값 유지" 전략
        // (txt에 없는 key는 건드리지 않음)

        std::string line;
        int updated = 0;

        auto set_if = [&](const std::string& key, double& var, double v){
            if (key == "") return;
            var = v;
            updated++;
        };

        while (std::getline(file, line))
        {
            line = trim_copy(line);
            if (line.empty() || line[0] == '#') continue;

            auto pos = line.find('=');
            if (pos == std::string::npos) continue;

            std::string key = trim_copy(line.substr(0, pos));
            std::string val = trim_copy(line.substr(pos + 1));

            try {
                double v = std::stod(val);

                if (key == "DP_x_rel") set_if(key, DP_x_rel, v);
                else if (key == "DP_y_rel") set_if(key, DP_y_rel, v);
                else if (key == "DP_psi_rel") set_if(key, DP_psi_rel, v);

                else if (key == "offset_x") set_if(key, offset_x_, v);
                else if (key == "offset_y") set_if(key, offset_y_, v);

                else if (key == "cradle_x_rel") set_if(key, cradle_x_rel, v);

                else if (key == "duration_time") set_if(key, duration_time, v);
                else if (key == "dp_dist") set_if(key, dp_dist, v);

                else if (key == "stream_speed") set_if(key, stream_speed, v);

                else if (key == "dock_radius") set_if(key, dock_radius, v);

                else if (key == "MS_speed") set_if(key, MS_speed, v);
                
                
                else if (key == "Xu") set_if(key, Xu, v);
                else if (key == "Xuu") set_if(key, Xuu, v);
                else if (key == "Yv") set_if(key, Yv, v);
                else if (key == "Yvv") set_if(key, Yvv, v);
                else if (key == "Yr") set_if(key, Yr, v);
                else if (key == "Yrr") set_if(key, Yrr, v);
                else if (key == "Yrrr") set_if(key, Yrrr, v);
                else if (key == "Nr") set_if(key, Nr, v);
                else if (key == "Nrr") set_if(key, Nrr, v);
                else if (key == "Nrrr") set_if(key, Nrrr, v);
                else if (key == "Nv") set_if(key, Nv, v);
                else if (key == "u_F") set_if(key, u_F, v);
                else if (key == "v_F") set_if(key, v_F, v);
                else if (key == "r_F") set_if(key, r_F, v);
                else if (key == "Xu_dot") set_if(key, Xu_dot, v);
                else if (key == "Yv_dot") set_if(key, Yv_dot, v);
                else if (key == "Nr_dot") set_if(key, Nr_dot, v);
                else if (key == "M") set_if(key, M, v);
                else if (key == "I") set_if(key, I, v);

                else if (key == "alpha_side") set_if(key, alpha_side, v);
                else if (key == "alpha_dock") set_if(key, alpha_dock, v);

                else if (key == "cbf_d") set_if(key, cbf_d, v);
                else if (key == "cbf_a") set_if(key, cbf_a, v);
                else if (key == "cbf_b") set_if(key, cbf_b, v);
                else if (key == "cbf_x0") set_if(key, cbf_x0, v);
                else if (key == "cbf_y0") set_if(key, cbf_y0, v);

                else if (key == "hp") set_if(key, hp, v);
                else if (key == "final_boost_radius") set_if(key, final_boost_radius, v);
                else if (key == "final_end_radius") set_if(key, final_end_radius, v);

                // ---- MPC Q ----
                else if (key == "Q_x")     set_if(key, Q_x, v);
                else if (key == "Q_y")     set_if(key, Q_y, v);
                else if (key == "Q_psi")   set_if(key, Q_psi, v);
                else if (key == "Q_u")     set_if(key, Q_u, v);
                else if (key == "Q_v")     set_if(key, Q_v, v);
                else if (key == "Q_r")     set_if(key, Q_r, v);
                else if (key == "Q_delta") set_if(key, Q_delta, v);
                else if (key == "Q_F")     set_if(key, Q_F, v);

                // ---- MPC QN ----
                else if (key == "QN_x")     set_if(key, QN_x, v);
                else if (key == "QN_y")     set_if(key, QN_y, v);
                else if (key == "QN_psi")   set_if(key, QN_psi, v);
                else if (key == "QN_u")     set_if(key, QN_u, v);
                else if (key == "QN_v")     set_if(key, QN_v, v);
                else if (key == "QN_r")     set_if(key, QN_r, v);
                else if (key == "QN_delta") set_if(key, QN_delta, v);
                else if (key == "QN_F")     set_if(key, QN_F, v);

                // ---- MPC R ----
                else if (key == "R_ddelta") set_if(key, R_ddelta, v);
                else if (key == "R_dF")     set_if(key, R_dF, v);


                else if (key == "gamma") set_if(key, gamma, v);
                else if (key == "lookahead")     set_if(key, lookahead, v);
                else if (key == "Kp_steer")     set_if(key, Kp_steer, v);
                else if (key == "Kd_steer")     set_if(key, Kd_steer, v);
                else if (key == "Ki_steer")     set_if(key, Ki_steer, v);
                else if (key == "Kp_F")     set_if(key, Kp_F, v);
                else if (key == "Ki_F")     set_if(key, Ki_F, v);


            
            }
            catch (...) {
                RCLCPP_WARN(this->get_logger(),
                    "[PARAM] Parse failed line: %s", line.c_str());
            }
        }

        RCLCPP_INFO(this->get_logger(),
            "[PARAM] Loaded txt. updated=%d (%s)", updated, path.c_str());

        return true;
    }


    bool loadBRSFromCSV(const std::string& path)
    {
        std::ifstream file(path);
        if (!file.is_open())
        {
            RCLCPP_ERROR(this->get_logger(),
                "[BRS] Cannot open CSV: %s", path.c_str());
            return false;
        }

        brs_.data.clear();
        brs_.data.reserve(
            brs_.Nx * brs_.Nx * brs_.Nx * brs_.Nx
        );

        std::string line;
        while (std::getline(file, line))
        {
            std::stringstream ss(line);
            std::string cell;
            while (std::getline(ss, cell, ','))
            {
                brs_.data.push_back(std::stod(cell));
            }
        }
        file.close();

        const size_t expected = std::pow(brs_.Nx, 4);
        if (brs_.data.size() != expected)
        {
            RCLCPP_FATAL(this->get_logger(),
                "[BRS] Size mismatch: %zu vs %zu",
                brs_.data.size(), expected);
            return false;
        }

        RCLCPP_INFO(this->get_logger(),
            "[BRS] Loaded %zu values", brs_.data.size());

        return true;
    }

    bool updateBRSExist()
    {
        // =====================================
        // 1) cradle 기준 상대좌표 (x, y only)
        // =====================================
        // double x_rel_abs = x_ + hp * cos(yaw_) - CD_x - cradle_x_rel*cos(CD_psi) - offset_x_;
        // double y_rel_abs = y_ + hp * sin(yaw_) - CD_y - cradle_x_rel*sin(CD_psi) - offset_y_;
    
        // double x_rel = x_rel_abs*cos(CD_psi) + y_rel_abs*sin(CD_psi);
        // double y_rel = x_rel_abs*sin(-CD_psi) + y_rel_abs*cos(CD_psi);
        // 1) 월드에서 “후크포인트(선수 hp)” 위치
        const double px = x_ + hp * std::cos(yaw_);
        const double py = y_ + hp * std::sin(yaw_);

        // 2) 크래들 기준점(월드) : CD 위치 + 크래들 오프셋(크래들 x축 방향)
        const double cx = (CD_x + offset_x_) + cradle_x_rel * std::cos(CD_psi);
        const double cy = (CD_y + offset_y_) + cradle_x_rel * std::sin(CD_psi);

        // 3) 월드 상대벡터
        const double dx = px - cx;
        const double dy = py - cy;

        // 4) 월드 -> CD frame (R(-CD_psi))
        const double c = std::cos(CD_psi);
        const double s = std::sin(CD_psi);

        const double x_rel =  c * dx + s * dy;
        const double y_rel = -s * dx + c * dy;

        // RCLCPP_INFO_THROTTLE(
        //     this->get_logger(),
        //     *this->get_clock(),
        //     2000,
        //     "[BRS] x_rel=%.2f y_rel=%.2f",
        //     x_rel, y_rel
        // );
    
        if ( sqrt(x_rel*x_rel + y_rel*y_rel) < dock_radius || x_rel > dock_radius )
        {
            // RCLCPP_INFO_THROTTLE(
            //     this->get_logger(),
            //     *this->get_clock(),
            //     2000,
            //     "[BRS] Recovery success"
            // );
            return true;
        }

        // =====================================
        // 2) BRS grid boundary check (x,y only)
        // =====================================
        if ( x_rel < brs_.xmin || x_rel > brs_.xmax ||
             y_rel < brs_.ymin || y_rel > brs_.ymax )
        {
            RCLCPP_INFO_THROTTLE(
                this->get_logger(),
                *this->get_clock(),
                2000,
                "[BRS] OUTSIDE BRS GRID"
            );
            return false;
        }
    
        // =====================================
        // 3) index 계산 (x,y only)
        // =====================================
        auto idx = [&](double x, double xmin, double dx) {
            return static_cast<int>(
                std::round((x - xmin) / dx)
            );
        };
    
        int ix = idx(x_rel, brs_.xmin, brs_.dx);
        int iy = idx(y_rel, brs_.ymin, brs_.dy);
    
        if (ix < 0 || ix >= brs_.Nx ||
            iy < 0 || iy >= brs_.Nx)
        {
            return false;
        }
    
        // =====================================
        // 4) BRS lookup
        //    min over (psi, w) dimensions
        // =====================================
        auto brsIndex = [&](int ix, int iy, int iz, int iw) {
            return (((iw * brs_.Nx + iz) * brs_.Nx + iy) * brs_.Nx + ix);
        };
    
        double phi_min = std::numeric_limits<double>::infinity();
    
        for (int iz = 0; iz < brs_.Nx; ++iz)        // psi dimension
        {
            for (int iw = 0; iw < brs_.Nx; ++iw)    // 4th dimension
            {
                phi_min = std::min(
                    phi_min,
                    brs_.data[ brsIndex(ix, iy, iz, iw) ]
                );
            }
        }
    
        // =====================================
        // 5) inside / outside BRS
        // =====================================
        return (phi_min < 0.0);
            
    }
    


    // ============================================================
    // EKF Callback
    // ============================================================
    void ekf_callback(const std_msgs::msg::Float64MultiArray::SharedPtr msg)
    {
        if (only_xypsi){
        // const double real_x_   = msg->data[0];
        // const double real_y_   = msg->data[1];
        // const double real_psi_ = msg->data[2];
        const double real_u_   = msg->data[3];
        const double real_v_   = msg->data[4];
        const double real_r_   = msg->data[5];


        if (msg->data.size() < 3) return;

        const double x_meas   = msg->data[0];
        const double y_meas   = msg->data[1];
        const double psi_meas = msg->data[2];
    
        auto now = this->get_clock()->now();
    
        if (!pose_time_init_)
        {
            pose_time_init_ = true;
            last_pose_time_ = now;

            x_   = x_meas;
            y_   = y_meas;
            yaw_ = psi_meas;
        
            // (선택) 초기 u,v,r를 0으로
            u_ = 0.0; v_ = 0.0; r_ = 0.0;
            
            updateStateVector();
    
    
            // RCLCPP_INFO(this->get_logger(),
            // "[MPC] u=%.2f v=%.2f r=%.2f real_u=%.2f real_v=%.2f real_r=%.2f dt = %.2f",
            // u_, v_, r_, real_u_, real_v_, real_r_ , dt);
    
        }
    
        const double dt = (now - last_pose_time_).seconds();
        last_pose_time_ = now;
    
        // dt가 비정상적으로 크거나 작으면 방어
        if (dt <= 1e-4 || dt > 2.0) return;
    
        // ---- (추천) EKF step ----
        uvr_ekf_.step(x_meas, y_meas, psi_meas, dt);
    
        // 네 MPC 상태에 반영
        x_   = msg->data[0];
        y_   = msg->data[1];
        yaw_ = yaw_unwrap(msg->data[2], yaw_);  // 필터된 psi
        u_   = uvr_ekf_.u();
        v_   = uvr_ekf_.v();
        r_   = uvr_ekf_.r();
    
        updateStateVector();


        // RCLCPP_INFO(this->get_logger(),
        // "[MPC] E_u=%.2f E_v=%.2f E_r=%.2f dt = %.2f",
        // u_ - real_u_, v_ - real_v_, r_ - real_r_ , dt);

        // RCLCPP_INFO(this->get_logger(),
        // "[MPC] u=%.2f v=%.2f r=%.2f real_u=%.2f real_v=%.2f real_r=%.2f dt = %.2f",
        // u_, v_, r_, real_u_, real_v_, real_r_ , dt);


        }
        else{
            if (msg->data.size() < 6) return;

            x_ = msg->data[0];
            y_ = msg->data[1];
            yaw_ = yaw_unwrap(msg->data[2], yaw_);
            u_ = msg->data[3];
            v_ = msg->data[4];
            r_ = msg->data[5];

            updateStateVector();
        }

    }    
    

    void lars_state_callback(const std_msgs::msg::Float64MultiArray::SharedPtr msg)
    {
        if (msg->data.size() < 7) return;

        CD_x = msg->data[0] - offset_x_;
        CD_y = msg->data[1] - offset_y_;
        CD_psi = msg->data[2];
        MS_x = msg->data[3] - offset_x_;
        MS_y = msg->data[4] - offset_y_;
        MS_psi = msg->data[5];
        // stream_speed = msg->data[6];
    }    

    void lars_mode_callback(const std_msgs::msg::Float64::SharedPtr msg)
    {
        lars_mode_ = msg->data;
    }    

    void runDOB()
    {
        // EKF state 기반 DOB 수행
        std::array<double,8> state_for_dob = {
            states_[0], states_[1], states_[2],
            states_[3], states_[4], states_[5],
            states_[6], states_[7]
        };

        DOB::updateDOB(
            state_for_dob,
            dob_state_,
            dob_dt_, 
            stream_speed,
            Xu, Xuu,
            Yv, Yvv, Yr, Yrr, Yrrr,
            Nr, Nrr, Nrrr, Nv, 
            u_F, v_F, r_F, 
            Xu_dot, Yv_dot, Nr_dot,
            M, I
        );

        
        // DOB publish (10 Hz)
        std_msgs::msg::Float64MultiArray dob_msg;
        dob_msg.data = {
            dob_state_.param_filtered[0],
            dob_state_.param_filtered[1],
            dob_state_.param_filtered[2]
        };

        dob_pub_->publish(dob_msg);
    }

    // ============================================================
    // MAIN MPC LOOP
    // ============================================================
    void run()
    {
        auto now = this->get_clock()->now();
        double elapsed = (now - start_time_).seconds();
        
        // loadParamsFromTxt(param_txt_path_);
        // 5초마다만 txt reload
        if ((now - last_param_update_time_).seconds() >= 5.0) {
            loadParamsFromTxt(param_txt_path_);
            applyCostWeights(); 
            last_param_update_time_ = now;
        }
        // =====================================================
        // 2) MAIN MPC LOOP (warm-up 이후 매 타임스텝 반복)
        // =====================================================

        // ---- Reachability 판단 ----
        BRS_exist_ = updateBRSExist();

        // ---- Mode check ----
        ModeCheck();
        checkWaypointHoldAndSwitch();

        RCLCPP_INFO_THROTTLE(
            this->get_logger(),
            *this->get_clock(),
            2000,   // ms (2초마다 출력)
            "[BRS CHECK] BRS_exist = %s",
            BRS_exist_ ? "TRUE" : "FALSE"
        );
        

        // double LOS = 0.0;
        // double gamma = 2.0;
        // double lookahead = 5.0;
        // double beta = 0.0;
        // double Kp = 0.0;
        // double Ki = 0.0;
        // double ye = 0.0;
        // double start_x = 0.0;
        // double start_y = 0.0;
        // double ye = 0.0;
    
        //// ALOS guidance ////

        // 1) 좌표 통일: 여기서는 "local(frame with offset removed)" 사용
        //    x_local, y_local은 states_[0], states_[1]와 동일한 좌표계
        const double x_local = x_ - offset_x_;
        const double y_local = y_ - offset_y_;

        // 2) 현재 목표 waypoint (이미 local 좌표계로 관리 중)
        const double x_t = wpt_x;
        const double y_t = wpt_y;

        // 3) 목표가 바뀌면(path segment가 바뀌면) start_x/start_y 재설정
        //    - docking처럼 목표가 움직이면 너무 자주 바뀔 수 있으니 threshold를 둠
        const double target_change = std::hypot(x_t - last_target_x_, y_t - last_target_y_);
        if (!alos_init_ || target_change > 0.5)  // 0.5m는 예시 (원하면 조절)
        {
            // 세그먼트 시작점을 "현재 위치"로 두는 방식(실전에서 안정적)
            start_x = x_local + 3;
            start_y = y_local;

            last_target_x_ = x_t;
            last_target_y_ = y_t;

            // 필요하면 beta를 리셋(급격한 목표 변경 시 peaking 방지)
            // beta = 0.0;

            alos_init_ = true;
        }

        // 4) 경로 방위각 pi_h (start -> target)
        const double seg_dx = x_t - start_x;
        const double seg_dy = y_t - start_y;
        const double seg_norm = std::hypot(seg_dx, seg_dy);

        // 세그먼트 길이가 너무 짧으면 계산 스킵

        const double pi_h = std::atan2(seg_dy, seg_dx);

        // 5) (x_e^p, y_e^p) = R^T(pi_h) * ([x;y] - [x_i;y_i])
        //    y_e^p가 cross-track error
        const double dx = x_local - start_x;
        const double dy = y_local - start_y;

        const double c = std::cos(pi_h);
        const double s = std::sin(pi_h);

        const double x_ep =  c * dx + s * dy;
        const double y_ep = -s * dx + c * dy;
    
        
        
        
        // 6) ALOS 적응법칙 (41): beta_hat_dot = gamma * y / sqrt(Delta^2 + y^2)
        const double Delta = std::max(1e-3, lookahead);  // lookahead는 양수
        const double denom = std::sqrt(Delta*Delta + y_ep*y_ep);
        const double beta_dot = gamma * (y_ep / denom);

        beta += beta_dot * dt_control_;

        // (선택) crab angle 추정치 제한 (작은 각도 가정 + 튜닝 안정성)
        beta = std::clamp(beta, -1.2, 1.2); // 약 +-34deg (원하면 조절)
        beta = -uvr_filter::wrapToPi(beta);

        // 7) ALOS 헤딩 지령 (40): psi_d = pi_h - beta_hat - atan(y/Delta)
        const double psi_d = uvr_filter::wrapToPi(pi_h - beta - std::atan2(y_ep, Delta) - yaw_);

        // 디버그(원하면)
        RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 500,
            "[ALOS] y_ep=%.2f beta_hat=%.3f psi_d=%.3f pi_h=%.3f",
            y_ep, beta, psi_d, pi_h);

        // RCLCPP_INFO_THROTTLE(this->get_logger(), *this->get_clock(), 500,
        // "[ALOS] start x=%.2f start y=%.3f last x=%.3f last y=%.3f",
        // start_x, start_y, last_target_x_, last_target_y_);
    


        const double dist = last_target_x_ - x_local ;

        dist_int += dist*dt_control_;

        if (dist_int < dist_int_max)
        {
            dist_int = dist_int_max;
        }
        else if (dist_int < -dist_int_max)
        {
            dist_int = -dist_int_max;
        }

        steer_int += psi_d*dt_control_;
        F_ = -Xuu*stream_speed*stream_speed/u_F + Kp_F * dist + Ki_F * dist_int;
        delta_ = -Kp_steer * psi_d - Kd_steer*r_ + Ki_steer*steer_int;

        //// ALOS guidance ////


        if (F_ >= 20.0){
            F_ = 20.0 - 0.1;
        }
        else if (F_ <= -5.0){
            F_ = -5 + 0.1;
        }
        
        if (delta_ > 6000*3.1415/180){
            delta_ = 6000*3.1415/180 - 0.2;
        }
        else if (delta_ < -6000*3.1415/180){
            delta_ = -6000*3.1415/180 + 0.2;
        }
        

        // Final_phase thrust
        if (final_get_in){
            double recovery_x = x_ - CD_x - offset_x_;
            double recovery_y = y_ - CD_y - offset_y_;

            if (!mission_over){
                if (sqrt(recovery_x*recovery_x + recovery_y*recovery_y) < final_end_radius){
                    mission_over = true;
                }
                else{
                    F_ = 20;
                    delta_ = 0.0;
                }
            }
            else{
                F_ = 0;
                delta_ = 0;
            }

        }

        // ---- PWM 변환 + Publish ----
        double delta_pwm = convertSteeringToPWM(delta_*0.01);
        double thrust_pwm = convertThrustToPWM(F_);
    
        std_msgs::msg::Float64MultiArray msg;
        msg.data = {delta_pwm, thrust_pwm, 0.0, 0.0};
        actuator_pub_->publish(msg);
    

        // ---- DEBUG PRINT ----
        // RCLCPP_INFO(this->get_logger(),
        //     "[MPC] x=%.2f y=%.2f yaw=%.2f u=%.2f v=%.2f r=%.2f  |  delta=%.2f F=%.2f",
        //     x_, y_, yaw_, u_, v_, r_, delta_, F_);

        // RCLCPP_INFO(this->get_logger(),
        // "[MPC] d_delta=%.2f d_throttle=%.2f", uMPC[0], uMPC[1]);
        
        // RCLCPP_INFO(this->get_logger(),
        // "[MPC state] x=%.2f y=%.2f, delta=%.2f F=%.2f", states_[0], states_[1], states_[6], states_[7]);      
        
        //         RCLCPP_INFO(this->get_logger(),
        // "[DOB] du=%.2f, dv=%.2f, dr=%.2f", dob_state_.param_filtered[0], dob_state_.param_filtered[1], dob_state_.param_filtered[2]);     
        // ---- DEBUG PRINT ----
    

        //--------------------------------------------
        // MPC VISUALIZATION MESSAGE (like Python)
        //--------------------------------------------
        aura_msg::msg::MPCTraj traj_msg;

        traj_msg.pred_num = N_;
        traj_msg.sampling_time = dt_control_;
        traj_msg.cpu_time = 0.0;   // acados stats 사용 원하면 수정 가능
        traj_msg.ref_num = 0.0;
        traj_msg.ref_dt = ref_dt_;

        traj_msg.traj_x = wpt_x;
        traj_msg.traj_y = wpt_y;
        traj_msg.theta = theta_;
        traj_msg.a = A_;
        traj_msg.b = B_;
        traj_msg.c = C_;

        // ----------------------------
        // predicted states (0~N)
        // ----------------------------
        for (int j = 0; j <= N_; j++)
        {
            double xj[HERON_NX];
            ocp_nlp_out_get(
                capsule_->nlp_config,
                capsule_->nlp_dims,
                capsule_->nlp_out,
                j,
                "x",
                xj);

            aura_msg::msg::MPCState st;
            st.x = xj[0] + offset_x_;
            st.y = xj[1] + offset_y_;
            st.p = xj[2];
            st.u = xj[3];
            st.v = xj[4];
            st.r = xj[5];
            st.delta = xj[6];
            st.f = xj[7];

            traj_msg.state.push_back(st);
        }

        // ----------------------------
        // reference (0~N)
        // ----------------------------
        for (int j = 0; j <= N_; j++)
        {
            aura_msg::msg::MPCState rf;
            rf.x = wpt_x;
            rf.y = wpt_y;
            rf.p = wpt_psi;
            rf.u = 0.0;
            rf.v = 0.0;
            rf.r = 0.0;
            rf.delta = 0.0;
            rf.f = 0.0;

            traj_msg.ref.push_back(rf);
        }

        // ----------------------------
        // obstacle info (for visualization)
        // ----------------------------
        double obs_vis[HERON_NP] = {
            1000.0, 40.0, 6.0,
            0.0,   1000.0, 8.0,
            0.0, 0.0, 0.0, 
            0.0, 0.0
        };

        for (int i = 0; i < 2; i++)
        {
            aura_msg::msg::ObsState ob;
            ob.x = obs_vis[i*3 + 0] + offset_x_;
            ob.y = obs_vis[i*3 + 1] + offset_y_;
            ob.rad = obs_vis[i*3 + 2];
            traj_msg.obs.push_back(ob);
        }


        mpcvis_pub_->publish(traj_msg);

    }
    


    // ============================================================
    // Reference Trajectory
    // ============================================================
    void generateReference()
    {
        reference_ = ReferenceGenerator::generateFigureEight(
            1000.0, ref_dt_, ref_x_, ref_y_, theta_, A_, B_, C_);

        std::vector<RefPoint> reduced;
        int step = int(dt_control_ / ref_dt_);

        for (int i = 0; i < reference_.size(); i += step)
            reduced.push_back(reference_[i]);

        reference_ = reduced;
    }


    // ============================================================
    // State helper functions
    // ============================================================
    void updateStateVector()
    {
        states_[0] = x_ - offset_x_;
        states_[1] = y_ - offset_y_;
        states_[2] = yaw_;
        states_[3] = u_;
        states_[4] = v_;
        states_[5] = r_;
        states_[6] = delta_;
        states_[7] = F_;

        hp_states_[0] = x_ - offset_x_;// + hp*cos(yaw_);
        hp_states_[1] = y_ - offset_y_;// + hp*sin(yaw_);
        hp_states_[2] = yaw_;
        hp_states_[3] = u_;
        hp_states_[4] = v_;
        hp_states_[5] = r_;
        hp_states_[6] = delta_;
        hp_states_[7] = F_;

        // RCLCPP_INFO(this->get_logger(),
        //     "[MPC] x=%.2f y=%.2f yaw=%.2f u=%.2f v=%.2f r=%.2f  |  delta=%.2f F=%.2f",
        //     states_[0], states_[1], states_[2], states_[3], states_[4], states_[5], states_[6], states_[7]);
    }

    void ModeCheck()
    {
        // =====================================
        // BRS 이탈 → 동조 기동 추종으로 복귀
        // =====================================
        if (!BRS_exist_ && lars_mode_ == 1)
        {
            lars_mode_ = 0.0;
            switch_1   = 0.0;
    
            lars_auto_triggered_ = false;
            inside_wpt_radius_  = false;
    
            // wpt_x   = CD_x + DP_x_rel;
            // wpt_y   = CD_y + DP_y_rel;
            // wpt_psi = CD_psi + DP_psi_rel;
            wpt_x   = CD_x + DP_x_rel*cos(CD_psi) - DP_y_rel*sin(CD_psi);
            wpt_y   = CD_y + DP_y_rel*cos(CD_psi) + DP_x_rel*sin(CD_psi);
            wpt_psi = CD_psi + DP_psi_rel;
            return;
        }
    
        // 정상 FSM
        if (lars_mode_ == 0.0)
        {
            // wpt_x   = CD_x + DP_x_rel;
            // wpt_y   = CD_y + DP_y_rel;
            // wpt_psi = CD_psi + DP_psi_rel;
            wpt_x   = CD_x + DP_x_rel*cos(CD_psi) - DP_y_rel*sin(CD_psi);
            wpt_y   = CD_y + DP_y_rel*cos(CD_psi) + DP_x_rel*sin(CD_psi);
            wpt_psi = CD_psi + DP_psi_rel;
            switch_1 = 0.0;
        }
        else if (lars_mode_ == 1.0)
        {
            wpt_x   = CD_x;
            wpt_y   = CD_y;
            wpt_psi = CD_psi;
            switch_1 = 1.0;
        }
        else if (lars_mode_ == 2.0)
        {
            wpt_x   = CD_x;
            wpt_y   = CD_y;
            wpt_psi = CD_psi;
            switch_1 = 1.0; // 또는 final 전용 switch
        }


        // RCLCPP_INFO(this->get_logger(),
        //     "wpt_x=%.2f wpt_y=%.2f",
        //     wpt_x, wpt_y);
    }
    
    void checkWaypointHoldAndSwitch()
    {
        /// Terminal acting ///
        double hp_x = x_ + hp*cos(yaw_) - offset_x_;
        double hp_y = y_ + hp*sin(yaw_) - offset_y_;
        // RCLCPP_WARN(this->get_logger(),
        // "dock_check : %.2f",
        // sqrt((hp_x - CD_x)*(hp_x - CD_x) + (hp_y - CD_y)*(hp_y - CD_y)));
        if (sqrt((hp_x - CD_x)*(hp_x - CD_x) + (hp_y - CD_y)*(hp_y - CD_y)) < final_boost_radius){
            lars_mode_ = 2;
            final_get_in = true;

            RCLCPP_WARN(this->get_logger(),
            "Lars_mode : %.2f",
            lars_mode_);
        }
        
        // already promoted → do nothing
        if (lars_auto_triggered_) return;
    
        // BRS 내부일 때만 전이 평가   
        if (!BRS_exist_) {
            inside_wpt_radius_ = false;
            return;
        }

        // distance in global frame
        double dx = x_ - (wpt_x + offset_x_);
        double dy = y_ - (wpt_y + offset_y_);
        double dist = std::sqrt(dx*dx + dy*dy);
    
        auto now = this->get_clock()->now();
    
        if ((dist <= dp_dist) && (lars_mode_ == 0))
        {
            if (!inside_wpt_radius_)
            {
                inside_wpt_radius_ = true;
                wpt_enter_time_ = now;
            }
            else
            {
                double duration = (now - wpt_enter_time_).seconds();
                if (duration >= duration_time)
                {
                    lars_mode_ = 1.0;
                    lars_auto_triggered_ = true;
    
                    RCLCPP_WARN(this->get_logger(),
                        "[AUTO MODE] Waypoint held for %.2f s → lars_mode = 1",
                        duration);
                }
            }
        }
        else
        {
            inside_wpt_radius_ = false;
        }


    }
    
    void applyCostWeights()
    {
        // ----------------------------------
        // Stage W (HERON_NY x HERON_NY)
        // y = [x,y,psi,u,v,r,delta,F, d_delta, d_F]
        // ----------------------------------
        std::vector<double> W(HERON_NY * HERON_NY, 0.0);
    
        auto set_diag = [&](int idx, double val){
            W[idx * HERON_NY + idx] = 2.0 * val; // Python의 2*diag와 동일
        };
    
        // Q (state 8) : index 0..7
        set_diag(0, Q_x);
        set_diag(1, Q_y);
        set_diag(2, Q_psi);
        set_diag(3, Q_u);
        set_diag(4, Q_v);
        set_diag(5, Q_r);
        set_diag(6, Q_delta);
        set_diag(7, Q_F);
    
        // R (input 2) : index 8..9
        set_diag(8, R_ddelta); // d_delta
        set_diag(9, R_dF);     // d_F
    
        for (int j = 0; j < N_; ++j)
        {
            ocp_nlp_cost_model_set(
                capsule_->nlp_config,
                capsule_->nlp_dims,
                capsule_->nlp_in,
                j,
                "W",
                W.data()
            );
        }
    
        // ----------------------------------
        // Terminal WN (HERON_NYN x HERON_NYN)
        // yN = [x,y,psi,u,v,r,delta,F]
        // ----------------------------------
        std::vector<double> WN(HERON_NYN * HERON_NYN, 0.0);
    
        auto set_diagN = [&](int idx, double val){
            WN[idx * HERON_NYN + idx] = 2.0 * val;
        };
    
        set_diagN(0, QN_x);
        set_diagN(1, QN_y);
        set_diagN(2, QN_psi);
        set_diagN(3, QN_u);
        set_diagN(4, QN_v);
        set_diagN(5, QN_r);
        set_diagN(6, QN_delta);
        set_diagN(7, QN_F);
    
        ocp_nlp_cost_model_set(
            capsule_->nlp_config,
            capsule_->nlp_dims,
            capsule_->nlp_in,
            N_,
            "W",
            WN.data()
        );
    }
    


    void setInitialStateBounds()
    {
        ocp_nlp_constraints_model_set(
            capsule_->nlp_config, capsule_->nlp_dims,
            capsule_->nlp_in, capsule_->nlp_out,
            0, "lbx", states_);

        ocp_nlp_constraints_model_set(
            capsule_->nlp_config, capsule_->nlp_dims,
            capsule_->nlp_in, capsule_->nlp_out,
            0, "ubx", states_);

        // RCLCPP_INFO(this->get_logger(),
        //     "delta=%.2f F=%.2f",
        //     states_[6], states_[7]);            

    }

    double yaw_unwrap(double yaw_new, double yaw_prev)
    {
        double diff = yaw_new - yaw_prev;
        if (diff > M_PI) yaw_new -= 2*M_PI;
        else if (diff < -M_PI) yaw_new += 2*M_PI;
        return yaw_new;
    }

    double convertSteeringToPWM(double steer)
    {

        return steer;
    }

    double convertThrustToPWM(double thrust)
    {
        
        double T = (thrust >= 0) ? 5000*std::sqrt(thrust/35) : -5000*std::sqrt(-thrust/35);
        
        return T;
    }

private:

    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr actuator_pub_;
    rclcpp::Subscription<std_msgs::msg::Float64MultiArray>::SharedPtr ekf_sub_;
    rclcpp::Subscription<std_msgs::msg::Float64MultiArray>::SharedPtr lars_state_sub_;
    rclcpp::Subscription<std_msgs::msg::Float64>::SharedPtr lars_mode_sub_;
    rclcpp::Publisher<aura_msg::msg::MPCTraj>::SharedPtr mpcvis_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr dob_pub_;
    rclcpp::TimerBase::SharedPtr timer_;
    rclcpp::TimerBase::SharedPtr dob_timer_;


    heron_solver_capsule *capsule_;

    // States
    double x_, y_, yaw_, u_, v_, r_;
    double ref_x_, ref_y_;
    double delta_, F_;
    double offset_x_, offset_y_;

    double dt_control_;
    int N_;

    std::vector<RefPoint> reference_;
    double A_, B_, C_, theta_;
    double ref_dt_;
    bool warmed_up_ = false;


    rclcpp::Time start_time_;

    double states_[HERON_NX];
    double hp_states_[HERON_NX];

    // DOB
    DOB::DOBState dob_state_;
    double dob_dt_ = 0.1;   // Python과 동일

    // Cradle, Mother ship
    double CD_x;
    double CD_y;
    double CD_psi;

    double MS_x;
    double MS_y;
    double MS_psi;
    double MS_speed;

    double lars_mode_;
    double switch_1;

    double stream_speed;

    double wpt_x;
    double wpt_y;
    double wpt_psi;

    double DP_x_rel;
    double DP_y_rel;
    double DP_psi_rel;

    double cradle_x_rel;
    double dock_radius;
    double duration_time;
    double dp_dist;

    double Xu;
    double Xuu;
    double Yv;
    double Yvv;
    double Yr;
    double Yrr;
    double Yrrr;
    double Nr;
    double Nrr;
    double Nrrr;
    double Nv;
    double u_F;
    double v_F;
    double r_F;

    double Xu_dot;
    double Yv_dot;
    double Nr_dot;
    double M;
    double I;

    double alpha_side;
    double alpha_dock;

    double cbf_d;
    double cbf_a;
    double cbf_b;
    double cbf_x0;
    double cbf_y0;

    bool only_xypsi;

    double Q_x=1e3, Q_y=1e3, Q_psi=1e5, Q_u=1e1, Q_v=1e3, Q_r=1e2, Q_delta=1e1, Q_F=1e-1;
    double QN_x=1e5, QN_y=1e5, QN_psi=1e4, QN_u=1e-2, QN_v=1e-1, QN_r=1e2, QN_delta=1e1, QN_F=1e-1;
    double R_ddelta=1e2, R_dF=1e-2;

    // ---- LOS ---- 
    double LOS = 0.0;
    double gamma = 0.3;
    double lookahead = 5.0;
    double beta = 0.0;
    double Kp = 0.0;
    double Ki = 0.0;
    double start_x = 0.0;
    double start_y = 0.0;
    double ye = 0.0;
    bool alos_init_ = false;
    double last_target_x_ = 0.0;   // 이전 waypoint (local)
    double last_target_y_ = 0.0;
    double dist_int = 0.0;
    double dist_int_max = 5.0;

    double Kp_steer = 50.0;
    double Kd_steer = 0.0; 
    double Ki_steer = 0.0; 
    double Kp_F = 0.1;
    double Ki_F = 0.00;
    double steer_int = 0.0;
    double steer_int_max = 5.0;

    std::string param_txt_path_ = "/home/user/aura_ws/src/lars/config/mpc_params.txt";

    rclcpp::Time last_param_update_time_;

    // ---- u,v,r 추정기(EKF 추천) ----
    uvr_filter::UVREKF uvr_ekf_;   // 또는 uvr_filter::UVRFilterLPF uvr_lpf_;

    // ---- dt 계산용 ----
    rclcpp::Time last_pose_time_;
    bool pose_time_init_ = false;

    // ---- Auto lars_mode promotion ----
    bool inside_wpt_radius_ = false;
    bool lars_auto_triggered_ = false;
    rclcpp::Time wpt_enter_time_;

    // Final mode
    bool final_get_in = false;
    bool mission_over = false;
    double final_boost_radius;
    double final_end_radius;

    double hp;

    // ===================================================
    // BRS GRID (Python 정의와 정합)
    // ===================================================
    struct BRSGrid {
        int Nx = 21;
        std::vector<double> data;   // size = Nx^4

        // Python mins/maxs와 동일
        double xmin = -10.0, xmax = 2.0;
        double ymin = -6.0, ymax = 6.0;
        double psimin = -3.141592, psimax = 3.141592;

        double dx, dy, dpsi;
    };

    BRSGrid brs_;
    bool BRS_exist_ = false;

};


int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<AuraMPC>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
