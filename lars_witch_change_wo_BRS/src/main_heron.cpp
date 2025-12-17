#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <std_msgs/msg/float64.hpp>

#include "reference_generator.hpp"
#include "dob.hpp"
#include "acados_solver_heron.h"

#include "aura_msg/msg/mpc_state.hpp"
#include "aura_msg/msg/mpc_traj.hpp"
#include "aura_msg/msg/obs_state.hpp"

#include <cmath>
#include <chrono>
#include <vector>

using namespace std::chrono_literals;

static constexpr double ship_state_x = 290245.0;
static constexpr double ship_state_y = 4118000.0;

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
        wpt_y = 0.0;
        wpt_psi = 0.0;

        DP_x_rel = -15.0;
        DP_y_rel = -5.0;
        DP_psi_rel = 0.0;


        N_ = HERON_N;
        dt_control_ = 0.5;

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

        timer_ = this->create_wall_timer(
            std::chrono::duration<double>(dt_control_),
            std::bind(&AuraMPC::run, this));

        dob_timer_ = this->create_wall_timer(
            std::chrono::duration<double>(0.1),   // 10 Hz
            std::bind(&AuraMPC::runDOB, this)
        );            

        start_time_ = this->get_clock()->now();

        RCLCPP_INFO(this->get_logger(), "AuraMPC initialized.");
    }

    ~AuraMPC()
    {
        heron_acados_free(capsule_);
        heron_acados_free_capsule(capsule_);
    }

private:

    // ============================================================
    // EKF Callback
    // ============================================================
    void ekf_callback(const std_msgs::msg::Float64MultiArray::SharedPtr msg)
    {
        if (msg->data.size() < 6) return;

        x_ = msg->data[0];
        y_ = msg->data[1];
        yaw_ = yaw_unwrap(msg->data[2], yaw_);
        u_ = msg->data[3];
        v_ = msg->data[4];
        r_ = msg->data[5];

        updateStateVector();
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
        stream_speed = msg->data[6];
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
            dob_dt_
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
    
        // ---------------------------------------
        // 1) WARM-UP (딱 1회만 실행)
        // ---------------------------------------
        if (!warmed_up_)
        {
            RCLCPP_INFO(this->get_logger(), "[WARM-UP] running warm-up cycle...");

            ref_x_ = x_; 
            ref_y_ = y_;

            ModeCheck();
            setInitialStateBounds(); // 초기 상태 고정
    
            int rti_phase = 1;
            ocp_nlp_solver_opts_set(
                capsule_->nlp_config,
                capsule_->nlp_opts,
                "rti_phase",
                &rti_phase);
    
            heron_acados_solve(capsule_);
    
            // warm-up 완료 표시
            warmed_up_ = true;
            RCLCPP_INFO(this->get_logger(), "[WARM-UP] completed.");
    
            return;  // warm-up 끝나고 다음 타임스텝부터 본 실행
        }
    
        // =====================================================
        // 2) MAIN MPC LOOP (warm-up 이후 매 타임스텝 반복)
        // =====================================================
        // ---- Mode check ----
        ModeCheck();
        checkWaypointHoldAndSwitch();

        // ---- 레퍼런스 설정 ----
        for (int j = 0; j <= N_; ++j)
        {
    
            double xr = wpt_x;
            double yr = wpt_y;
            double pr = wpt_psi;
    
            if (j == N_) {
                double yrefN[HERON_NYN] = {xr, yr, pr, 0.0, 0,0,0,0};
                ocp_nlp_cost_model_set(
                    capsule_->nlp_config,
                    capsule_->nlp_dims,
                    capsule_->nlp_in,
                    j,
                    "yref",
                    yrefN);
            } else {
                double yref[HERON_NY] = {xr, yr, pr, 0.0, 0,0,0,0,0,0};
                ocp_nlp_cost_model_set(
                    capsule_->nlp_config,
                    capsule_->nlp_dims,
                    capsule_->nlp_in,
                    j,
                    "yref",
                    yref);
            }

            // RCLCPP_INFO(this->get_logger(),
            // "[MPC state] xr=%.2f yr=%.2f", xr, yr);     
            
        }
    
        
        // ---- 장애물 파라미터 설정 ----
        double obs[HERON_NP] = {
            CD_x, CD_y, CD_psi,
            MS_x, MS_y, MS_psi,
            dob_state_.param_filtered[0],  dob_state_.param_filtered[1],  dob_state_.param_filtered[2],
            stream_speed, switch_1
        };
    
        for (int j = 0; j <= N_; j++)
            heron_acados_update_params(capsule_, j, obs, HERON_NP);
    
        // ---- RTI Preparation ----
        int rti_phase = 1;
        ocp_nlp_solver_opts_set(
            capsule_->nlp_config,
            capsule_->nlp_opts,
            "rti_phase",
            &rti_phase);
    
        heron_acados_solve(capsule_);
    
        // ---- state equality constraint update ----
        setInitialStateBounds();
    
        // ---- RTI Feedback ----
        rti_phase = 2;
        ocp_nlp_solver_opts_set(
            capsule_->nlp_config,
            capsule_->nlp_opts,
            "rti_phase",
            &rti_phase);
    
        heron_acados_solve(capsule_);
    


        // ---- Control output ----
        double uMPC[HERON_NU];
        ocp_nlp_out_get(capsule_->nlp_config, capsule_->nlp_dims,
                        capsule_->nlp_out, 0, "u", uMPC);
    
        delta_ += uMPC[0] * dt_control_;
        F_     += uMPC[1] * dt_control_;

        if (F_ > 40.0){
            F_ = 40.0;
        }

        
        // ---- PWM 변환 + Publish ----
        double delta_pwm = convertSteeringToPWM(delta_);
        double thrust_pwm = convertThrustToPWM(F_ * 100.0);
    
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
    }

    void ModeCheck()
    {
        if (lars_mode_ == 0.0){
            wpt_x = CD_x + DP_x_rel;
            wpt_y = CD_y + DP_y_rel;
            wpt_psi = CD_psi + DP_psi_rel;

            switch_1 = 0.0; 
        }
        else if (lars_mode_ == 1.0){
            wpt_x = CD_x;
            wpt_y = CD_y;
            wpt_psi = CD_psi;

            switch_1 = 1.0; 
        }
    }    
    void checkWaypointHoldAndSwitch()
    {
        // already promoted → do nothing
        if (lars_auto_triggered_) return;
    
        // distance in global frame
        double dx = x_ - (wpt_x + offset_x_);
        double dy = y_ - (wpt_y + offset_y_);
        double dist = std::sqrt(dx*dx + dy*dy);
    
        auto now = this->get_clock()->now();
    
        if (dist <= 1.5)
        {
            if (!inside_wpt_radius_)
            {
                inside_wpt_radius_ = true;
                wpt_enter_time_ = now;
            }
            else
            {
                double duration = (now - wpt_enter_time_).seconds();
                if (duration >= 5.0)
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
        //     "[MPC] x=%.2f y=%.2f yaw=%.2f u=%.2f v=%.2f r=%.2f  |  delta=%.2f F=%.2f",
        //     states_[0], states_[1], states_[2], states_[3], states_[4], states_[5], states_[6], states_[7]);            

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
        if (steer >= 300.0) return 2000.0;
        if (steer <= -300.0) return 1000.0;
        if (steer > 0.0 && steer < 300.0) return 1550.0 + steer * 1.6667;
        if (steer > -300.0 && steer < 0.0) return 1450.0 + steer * 1.6667;

        return 1500.0 + steer * 1.6667;
    }

    double convertThrustToPWM(double thrust)
    {

        double T = (thrust >= 0) ? std::sqrt(thrust) : -std::sqrt(-thrust);

        if (T >= 0)
            return std::clamp(1550.0 + 3.9*T, 1550.0, 2000.0);
        else
            return std::clamp(1450.0 + 3.9*T, 1000.0, 1450.0);
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

    double lars_mode_;
    double switch_1;

    double stream_speed;

    double wpt_x;
    double wpt_y;
    double wpt_psi;

    double DP_x_rel;
    double DP_y_rel;
    double DP_psi_rel;

    // ---- Auto lars_mode promotion ----
    bool inside_wpt_radius_ = false;
    bool lars_auto_triggered_ = false;
    rclcpp::Time wpt_enter_time_;
};


int main(int argc, char *argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<AuraMPC>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
