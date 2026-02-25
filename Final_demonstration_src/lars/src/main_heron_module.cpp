#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <std_msgs/msg/float64.hpp>

#include "uvr_filter.hpp"   // now includes PoseLPFEKF
#include "dob.hpp"          // now includes DobParams + loadDobParamsFromTxt
#include "mpc.hpp"          // new MPC module

#include "aura_msg/msg/mpc_state.hpp"
#include "aura_msg/msg/mpc_traj.hpp"
#include "aura_msg/msg/obs_state.hpp"

#include <rclcpp/rclcpp.hpp>
#include <rmw/qos_profiles.h>   // (필수는 아니지만 명확히)
// ============================================================
// ROS2 Wrapper Node
// - Reads state from ROS2 topics
// - Runs LPF+EKF -> (x,y,psi,u,v,r)
// - Runs DOB (10 Hz)
// - Runs MPC (dt_control)
// - Publishes actuator command
// ============================================================

class AuraMPCNode : public rclcpp::Node
{
public:
    AuraMPCNode()
    : Node("aura_mpc")
    {
        // -----------------------------
        // ROS2 I/O
        // -----------------------------
        auto qos_best_effort = rclcpp::QoS(rclcpp::KeepLast(10))
        .reliability(RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT)
        .durability(RMW_QOS_POLICY_DURABILITY_VOLATILE);

        
        actuator_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>(
            "/USVControlType", qos_best_effort);

        dob_pub_ = this->create_publisher<std_msgs::msg::Float64MultiArray>(
            "/DOB", qos_best_effort);

        ekf_sub_ = this->create_subscription<std_msgs::msg::Float64MultiArray>(
            "/StatusType", qos_best_effort,
            std::bind(&AuraMPCNode::ekf_callback, this, std::placeholders::_1));

        lars_state_sub_ = this->create_subscription<std_msgs::msg::Float64MultiArray>(
            "/lars_state", qos_best_effort,
            std::bind(&AuraMPCNode::lars_state_callback, this, std::placeholders::_1));

        mpcvis_pub_ = this->create_publisher<aura_msg::msg::MPCTraj>("/mpc_vis", qos_best_effort);
        // -----------------------------
        // Parameters
        // -----------------------------
        // NOTE: keep same path as your original. Change if needed.
        param_txt_path_ = "/home/user/aura_ws/src/lars/config/mpc_params_8.txt";
        brs_csv_path_   = "/home/user/aura_ws/src/lars/BRS_file/brs.csv";

        // Load DOB params (stream + model)
        if (!DOB::loadDobParamsFromTxt(param_txt_path_, dob_params_))
        {
            RCLCPP_WARN(this->get_logger(), "[DOB] Failed to load txt: %s", param_txt_path_.c_str());
        }

        // Create MPC module
        if (!mpc_.init(param_txt_path_, brs_csv_path_))
        {
            RCLCPP_ERROR(this->get_logger(), "[MPC] init failed");
        }

        // -----------------------------
        // Timers
        // -----------------------------
        // DOB @10Hz
        dob_timer_ = this->create_wall_timer(
            std::chrono::duration<double>(dob_dt_),
            std::bind(&AuraMPCNode::runDOB, this));

        // MPC @ dt_control (from module)
        mpc_timer_ = this->create_wall_timer(
            std::chrono::duration<double>(mpc_.dt_control()),
            std::bind(&AuraMPCNode::runMPC, this));

        RCLCPP_INFO(this->get_logger(), "AuraMPCNode started. dt_control=%.3f dob_dt=%.3f", mpc_.dt_control(), dob_dt_);
    }

private:
    // ============================================================
    // Callbacks
    // ============================================================
    void ekf_callback(const std_msgs::msg::Float64MultiArray::SharedPtr msg)
    {
        if (msg->data.size() < 3) return;

        const double x_meas = msg->data[0];
        const double y_meas = msg->data[1];
        const double psi_meas = msg->data[2];

        const double t = this->get_clock()->now().seconds();

        uvr_filter::PoseLPFEKF::Output out;
        if (!pose_est_.updatePose(x_meas, y_meas, psi_meas, t, out))
            return;

        // Save latest filtered/estimated state
        est_x_ = out.x;
        est_y_ = out.y;
        est_psi_ = out.psi;
        est_u_ = out.u;
        est_v_ = out.v;
        est_r_ = out.r;
        have_state_ = true;

        // Feed MPC module (it keeps internal delta,F)
        mpc_.setEstimatedState(est_x_, est_y_, est_psi_, est_u_, est_v_, est_r_);
    }

    void lars_state_callback(const std_msgs::msg::Float64MultiArray::SharedPtr msg)
    {
        if (msg->data.size() < 6) return;

        MPCModule::LarsState ls;
        ls.CD_x   = msg->data[0];
        ls.CD_y   = msg->data[1];
        ls.CD_psi = msg->data[2];
        ls.MS_x   = msg->data[3];
        ls.MS_y   = msg->data[4];
        ls.MS_psi = msg->data[5];
        mpc_.setLarsState(ls);
    }


    // ============================================================
    // DOB loop (10 Hz)
    // ============================================================
    void runDOB()
    {
        if (!have_state_) return;

        const double now_sec = this->get_clock()->now().seconds();

        // reload DOB txt every 5 seconds
        if (!dob_time_init_)
        {
            dob_time_init_ = true;
            last_dob_reload_sec_ = now_sec;
        }
        if ((now_sec - last_dob_reload_sec_) >= 5.0)
        {
            (void)DOB::loadDobParamsFromTxt(param_txt_path_, dob_params_);
            last_dob_reload_sec_ = now_sec;
        }

        // Build input for DOB
        std::array<double,8> state_for_dob = {
            est_x_, est_y_, est_psi_,
            est_u_, est_v_, est_r_,
            mpc_.delta(), mpc_.F()
        };

        DOB::updateDOB(
            state_for_dob,
            dob_state_,
            dob_dt_,
            dob_params_.stream_speed,
            dob_params_.Xu, dob_params_.Xuu,
            dob_params_.Yv, dob_params_.Yvv, dob_params_.Yr, dob_params_.Yrr, dob_params_.Yrrr,
            dob_params_.Nr, dob_params_.Nrr, dob_params_.Nrrr, dob_params_.Nv,
            dob_params_.u_F, dob_params_.v_F, dob_params_.r_F,
            dob_params_.Xu_dot, dob_params_.Yv_dot, dob_params_.Nr_dot,
            dob_params_.M, dob_params_.I
        );

        // expose filtered disturbances
        dob_u_ = dob_state_.param_filtered[0];
        dob_v_ = dob_state_.param_filtered[1];
        dob_r_ = dob_state_.param_filtered[2];

        // Publish DOB for debugging
        std_msgs::msg::Float64MultiArray msg;
        msg.data = {dob_u_, dob_v_, dob_r_};
        dob_pub_->publish(msg);
    }

    // ============================================================
    // MPC loop (dt_control)
    // ============================================================
    void runMPC()
    {
        if (!have_state_) return;

        const double now_sec = this->get_clock()->now().seconds();

        double new_delta = mpc_.delta();
        double new_F     = mpc_.F();

        // Use DOB outputs (already filtered)
        const bool ok = mpc_.step(dob_u_, dob_v_, dob_r_, now_sec, new_delta, new_F);
        if (!ok) return;

        // std::fprintf(stderr,
        //     "[Control] F=%.3f delta=%.3f\n",
        //     new_F, new_delta);
        
        mpc_.setActuatorInternal(new_delta, new_F);

        // Convert to PWM (keep your original conversion)
        const double delta_pwm  = convertSteeringToPWM(new_delta * 0.01);
        const double thrust_pwm = convertThrustToPWM(new_F);

        std_msgs::msg::Float64MultiArray cmd;
        cmd.data = {thrust_pwm, -delta_pwm, 0.0};
        actuator_pub_->publish(cmd);


        aura_msg::msg::MPCTraj traj_msg;
        traj_msg.pred_num = mpc_.N();
        traj_msg.sampling_time = mpc_.dt_control();
        traj_msg.cpu_time = 0.0;
        traj_msg.ref_num = 0.0;
        traj_msg.ref_dt = 0.0;          // main에서 관리 중인 값 사용

        traj_msg.traj_x = 0.0;            // main에서 관리 중인 값 사용
        traj_msg.traj_y = 0.0;
        traj_msg.theta = 0.0;
        traj_msg.a = 100.0;
        traj_msg.b = 100.0;
        traj_msg.c = 100.0;
        
        // ----------------------------
        // predicted states (0~N) : x,y만
        // ----------------------------
        std::vector<double> pred_x, pred_y;
        if (mpc_.getPredictedXY(pred_x, pred_y, /*absolute=*/true))
        {
            traj_msg.state.reserve(pred_x.size());
            for (size_t j = 0; j < pred_x.size(); ++j)
            {
                aura_msg::msg::MPCState st;
                st.x = pred_x[j];
                st.y = pred_y[j];
                traj_msg.state.push_back(st);
            }
        }

        // publish
        mpcvis_pub_->publish(traj_msg);



    }

    // ============================================================
    // PWM conversion (keep original behavior)
    // ============================================================
    static double convertSteeringToPWM(double steer)
    {
        // return steer;
        return (10000.0 * 180.0 / (30.0 * 3.141592)) * steer;
    }

    static double convertThrustToPWM(double thrust)
    {
        const double T = (thrust >= 0.0)
            ? 1000.0 * std::sqrt(thrust)
            : -1000.0 * std::sqrt(-thrust);
        return T;
    }

private:
    // ROS pubs/subs
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr actuator_pub_;
    rclcpp::Publisher<std_msgs::msg::Float64MultiArray>::SharedPtr dob_pub_;
    rclcpp::Publisher<aura_msg::msg::MPCTraj>::SharedPtr mpcvis_pub_;

    rclcpp::Subscription<std_msgs::msg::Float64MultiArray>::SharedPtr ekf_sub_;
    rclcpp::Subscription<std_msgs::msg::Float64MultiArray>::SharedPtr lars_state_sub_;

    rclcpp::TimerBase::SharedPtr dob_timer_;
    rclcpp::TimerBase::SharedPtr mpc_timer_;

    // Estimator (LPF + EKF)
    uvr_filter::PoseLPFEKF pose_est_;

    // Latest estimated state
    bool have_state_ = false;
    double est_x_ = 0.0, est_y_ = 0.0, est_psi_ = 0.0;
    double est_u_ = 0.0, est_v_ = 0.0, est_r_ = 0.0;

    // DOB
    DOB::DobParams dob_params_;
    DOB::DOBState dob_state_;
    double dob_dt_ = 0.1;

    double dob_u_ = 0.0, dob_v_ = 0.0, dob_r_ = 0.0;

    bool dob_time_init_ = false;
    double last_dob_reload_sec_ = 0.0;

    // MPC
    MPCModule mpc_;

    // paths
    std::string param_txt_path_;
    std::string brs_csv_path_;
};

int main(int argc, char* argv[])
{
    rclcpp::init(argc, argv);
    auto node = std::make_shared<AuraMPCNode>();
    rclcpp::spin(node);
    rclcpp::shutdown();
    return 0;
}
