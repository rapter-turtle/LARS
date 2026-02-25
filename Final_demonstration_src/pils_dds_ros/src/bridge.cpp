#include "pils_dds_bridge/bridge.hpp"
#include <chrono>

namespace pils_dds_bridge
{
BridgeNode::BridgeNode(const rclcpp::NodeOptions& options)
: rclcpp::Node("pils_dds_bridge_node", options),
  participant_(0),
  dds_sub_(participant_),
  dds_pub_(participant_),
  dds_topic_control_(dds::core::null),
  dds_writer_control_(dds::core::null),
  dds_topic_gwp_(dds::core::null),
  dds_reader_gwp_(dds::core::null),
  dds_topic_status_(dds::core::null),
  dds_reader_status_(dds::core::null),
  dds_reader_control_(dds::core::null),
  dds_writer_status_(dds::core::null), 
  dds_topic_cmd_init_(dds::core::null),
  dds_reader_cmd_init_(dds::core::null)
{ 
  init_parameters();

  participant_ = dds::domain::DomainParticipant(domain_id_);
  dds_sub_ = dds::sub::Subscriber(participant_);
  dds_pub_ = dds::pub::Publisher(participant_);

  init_dds();
  init_ros();

  const auto period = std::chrono::duration_cast<std::chrono::nanoseconds>(
    std::chrono::duration<double>(1.0 / rate_hz_));
  timer_ = create_wall_timer(period, std::bind(&BridgeNode::tick_10hz, this));

  RCLCPP_INFO(get_logger(),
              "독립 DDS<->ROS2 브리지 시작: domain_id=%d, rate=%.2fHz",
              domain_id_, rate_hz_);
}

void BridgeNode::init_parameters()
{
  declare_parameter<int>("domain_id", 200);
  declare_parameter<double>("rate_hz", 10.0);

  domain_id_ = get_parameter("domain_id").as_int();
  rate_hz_   = get_parameter("rate_hz").as_double();
}

void BridgeNode::init_ros()
{
  // auto qos = rclcpp::QoS(rclcpp::KeepLast(10));
  auto qos = rclcpp::QoS(rclcpp::KeepLast(10)).best_effort();  

  ros_sub_control_ = create_subscription<RosArray>(
    topic_usv_control_, qos,
    std::bind(&BridgeNode::on_ros_control, this, std::placeholders::_1));

  ros_pub_gwp_ = create_publisher<RosArray>(topic_global_waypoint_, qos);
  ros_pub_status_ = create_publisher<RosArray>(topic_status_, qos);

  // ros_sub_control_ = create_subscription<RosArray>(
  //   "USVControlType_ros", qos,
  //   std::bind(&BridgeNode::on_ros_control, this, std::placeholders::_1));
  // ros_pub_gwp_ = create_publisher<RosArray>("GlobalWaypointType_ros", qos);
  // ros_pub_status_ = create_publisher<RosArray>("StatusType_ros", qos);


    // ✅ 추가: ROS->DDS status
  ros_sub_status_ = create_subscription<RosArray>(
    topic_status_, qos,
    std::bind(&BridgeNode::on_ros_status, this, std::placeholders::_1));
  // ✅ 추가: DDS->ROS control
  ros_pub_control_ = create_publisher<RosArray>(topic_usv_control_, qos);

  ros_pub_cmd_init_ = create_publisher<RosArray>(topic_command_initial_, qos);

}

void BridgeNode::init_dds()
{
  // ROS 토픽 이름과 동일한 DDS topic name으로 구성
  dds_topic_control_ = dds::topic::Topic<DdsControl>(participant_, topic_usv_control_);
  dds_writer_control_ = dds::pub::DataWriter<DdsControl>(dds_pub_, dds_topic_control_);

  dds_topic_gwp_ = dds::topic::Topic<DdsGwp>(participant_, topic_global_waypoint_);
  dds_reader_gwp_ = dds::sub::DataReader<DdsGwp>(dds_sub_, dds_topic_gwp_);

  dds_topic_status_ = dds::topic::Topic<DdsStatus>(participant_, topic_status_);
  dds_reader_status_ = dds::sub::DataReader<DdsStatus>(dds_sub_, dds_topic_status_);

  // ✅ 추가: DDS->ROS control reader
  dds_reader_control_ = dds::sub::DataReader<DdsControl>(dds_sub_, dds_topic_control_);

  // ✅ 추가: ROS->DDS status writer
  dds_writer_status_ = dds::pub::DataWriter<DdsStatus>(dds_pub_, dds_topic_status_);

  dds_topic_cmd_init_ = dds::topic::Topic<DdsCmdInit>(participant_, topic_command_initial_);
  dds_reader_cmd_init_ = dds::sub::DataReader<DdsCmdInit>(dds_sub_, dds_topic_cmd_init_);


}

void BridgeNode::on_ros_control(const RosArray::SharedPtr msg)
{
  std::lock_guard<std::mutex> lock(m_ros_);
  last_control_ros_ = *msg;

  control_dirty_ = true;
  last_control_time_ = this->now();

  // publish_ros_control_to_dds();
}

void BridgeNode::tick_10hz()
{
  publish_ros_control_to_dds();
  poll_dds_and_publish_status();

  poll_dds_and_publish_cmd_init();

  poll_dds_and_publish_gwp();


  // publish_ros_status_to_dds();     
  // poll_dds_and_publish_control();  

}

// void BridgeNode::publish_ros_control_to_dds()
// {
//   std::optional<RosArray> copy;
//   {
//     std::lock_guard<std::mutex> lock(m_ros_);
//     copy = last_control_ros_;
//   }
//   if (!copy.has_value()) return;

//   const DdsControl dds_msg = ros_to_dds_control(*copy);
//   dds_writer_control_.write(dds_msg);

//   RCLCPP_INFO(get_logger(),
//   "[ROS->DDS] USVControlType write: rpm=%d steer=%d bow=%d",
//   (int)dds_msg.integratedTargetRPM(),
//   (int)dds_msg.integratedTargetSteering(),
//   (int)dds_msg.integratedBowThrust());
// }
void BridgeNode::publish_ros_control_to_dds()
{
  std::optional<RosArray> copy;
  rclcpp::Time stamp;
  {
    std::lock_guard<std::mutex> lock(m_ros_);
    if (!last_control_ros_.has_value()) return;
    if (!control_dirty_) return;          // ★ 핵심: 새 메시지 없으면 보내지 않음
    copy = last_control_ros_;
    stamp = last_control_time_;
    control_dirty_ = false;               // ★ 한번 보내면 dirty 해제
  }

  // 데이터 길이 체크 (비어있거나 2개 미만이면 0으로 매핑되어 나가버림)
  if (copy->data.size() < 2) {
    RCLCPP_WARN(get_logger(), "[ROS->DDS] USVControlType ignored: data.size()=%zu", copy->data.size());
    return;
  }

  const DdsControl dds_msg = ros_to_dds_control(*copy);
  dds_writer_control_.write(dds_msg);

  RCLCPP_INFO(get_logger(),
    "[ROS->DDS] USVControlType write: rpm=%d steer=%d bow=%d",
    (int)dds_msg.integratedTargetRPM(),
    (int)dds_msg.integratedTargetSteering(),
    (int)dds_msg.integratedBowThrust());
}


void BridgeNode::poll_dds_and_publish_gwp()
{
  DdsGwp latest;
  bool has = false;

  auto samples = dds_reader_gwp_.take();
  for (const auto& s : samples)
  {
    if (!s.info().valid()) continue;
    latest = s.data();
    has = true;
  }
  if (!has) return;

  RosArray ros = dds_to_ros_gwp(latest);
  ros_pub_gwp_->publish(ros);

}

void BridgeNode::poll_dds_and_publish_status()
{
  DdsStatus latest;
  bool has = false;

  auto samples = dds_reader_status_.take();
  for (const auto& s : samples)
  {
    if (!s.info().valid()) continue;
    latest = s.data();
    has = true;
  }
  if (!has) return;

  RosArray ros = dds_to_ros_status(latest);
  ros_pub_status_->publish(ros);

  RCLCPP_INFO(get_logger(),
  "[DDS->ROS] StatusType publish: size=%zu USV(x,y,h, u, v, r)=(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)",
  ros.data.size(),
  ros.data[0], ros.data[1], ros.data[2], ros.data[3], ros.data[4], ros.data[5]);
}


// ---------------------------
// actuator_outputs (ROS) -> USVControlType (DDS)
// ROS array: [steer_pwm, thrust_pwm, bow, ...]
// DDS fields:
//  - integratedTargetSteering <= steer_pwm
//  - integratedTargetRPM      <= thrust_pwm
//  - integratedBowThrust      <= bow
// ---------------------------
inline DdsControl ros_actuator_outputs_to_dds_control(const RosArray& ros)
{
  DdsControl dds;
  const auto& a = ros.data;

  const double steer = (a.size() > 0) ? a[0] : 0.0;
  const double thrust = (a.size() > 1) ? a[1] : 0.0;
  const double bow = (a.size() > 2) ? a[2] : 0.0;

  dds.integratedTargetSteering(clamp_to_short(steer));
  dds.integratedTargetRPM(clamp_to_short(thrust));
  dds.integratedBowThrust(clamp_to_short(bow));
  return dds;
}

// ---------------------------
// StatusType (DDS) -> /ekf/estimated_state (ROS)
// ROS array: [x, y, psi, u, v, r]  (size=6)
// u,v,r are estimated from (x,y,psi) using PoseLPFEKF (same as dds_to_ros_status)
// ---------------------------
inline RosArray dds_to_ros_ekf_estimated_state(const DdsStatus& s)
{
  RosArray ros;
  ros.data.resize(6);

  const double x = s.USV_x();
  const double y = s.USV_y();
  const double psi = s.USV_h(); // heading

  // time accumulation (same pattern as your dds_to_ros_status)
  using clock = std::chrono::steady_clock;
  static bool inited = false;
  static clock::time_point last_tp;
  static double t_sec = 0.0;

  const auto now_tp = clock::now();
  double dt = 0.0;

  if (!inited) {
    inited = true;
    last_tp = now_tp;
    dt = 0.0;
  } else {
    dt = std::chrono::duration<double>(now_tp - last_tp).count();
    last_tp = now_tp;
  }
  if (dt < 0.0) dt = 0.0;
  if (dt > 0.5) dt = 0.5;
  t_sec += dt;

  static uvr_filter::PoseLPFEKF pose_est_;
  uvr_filter::PoseLPFEKF::Output out;

  // 실패하면 (x,y,psi)만이라도 내보내고 u,v,r=0
  if (!pose_est_.updatePose(x, y, uvr_filter::wrapToPi(psi), t_sec, out)) {
    ros.data[0] = x;
    ros.data[1] = y;
    ros.data[2] = psi;
    ros.data[3] = 0.0;
    ros.data[4] = 0.0;
    ros.data[5] = 0.0;
    return ros;
  }

  ros.data[0] = x;
  ros.data[1] = y;
  ros.data[2] = psi;
  ros.data[3] = out.u;
  ros.data[4] = out.v; // keep your existing sign convention
  ros.data[5] = out.r;

  return ros;
}

static bool same_array(const std_msgs::msg::Float64MultiArray& a,
                       const std_msgs::msg::Float64MultiArray& b,
                       double eps = 1e-12)
{
  if (a.data.size() != b.data.size()) return false;
  for (std::size_t i = 0; i < a.data.size(); ++i) {
    if (std::fabs(a.data[i] - b.data[i]) > eps) return false;
  }
  return true;
}

void BridgeNode::on_ros_status(const RosArray::SharedPtr msg)
{
  std::lock_guard<std::mutex> lock(m_ros_);

  // ✅ echo loop 방지: DDS->ROS로 방금 뿌린 StatusType면 무시
  if (last_status_from_dds_ros_.has_value() && same_array(*msg, *last_status_from_dds_ros_)) {
    return;
  }

  last_status_ros_ = *msg;
  status_dirty_ = true;
  last_status_time_ = this->now();
}

void BridgeNode::publish_ros_status_to_dds()
{
  std::optional<RosArray> copy;
  {
    std::lock_guard<std::mutex> lock(m_ros_);
    if (!last_status_ros_.has_value()) return;
    if (!status_dirty_) return;
    copy = last_status_ros_;
    status_dirty_ = false;
  }

  // 최소 크기 체크 (StatusType은 보통 18 기대)
  if (copy->data.size() < 18) {
    RCLCPP_WARN(get_logger(), "[ROS->DDS] StatusType ignored: data.size()=%zu", copy->data.size());
    return;
  }

  const DdsStatus dds_msg = ros_to_dds_status(*copy);
  dds_writer_status_.write(dds_msg);

  RCLCPP_INFO(get_logger(),
    "[ROS->DDS] StatusType write: USV(x,y,h)=(%.3f, %.3f, %.3f)",
    copy->data[0], copy->data[1], copy->data[2]);
}







// static bool same_array(const std_msgs::msg::Float64MultiArray& a,
//                        const std_msgs::msg::Float64MultiArray& b,
//                        double eps = 1e-12)
// {
//   if (a.data.size() != b.data.size()) return false;
//   for (std::size_t i = 0; i < a.data.size(); ++i) {
//     if (std::fabs(a.data[i] - b.data[i]) > eps) return false;
//   }
//   return true;
// }

// void BridgeNode::on_ros_status(const RosArray::SharedPtr msg)
// {
//   std::lock_guard<std::mutex> lock(m_ros_);

//   // ✅ echo loop 방지: DDS->ROS로 방금 뿌린 StatusType면 무시
//   if (last_status_from_dds_ros_.has_value() && same_array(*msg, *last_status_from_dds_ros_)) {
//     return;
//   }

//   last_status_ros_ = *msg;
//   status_dirty_ = true;
//   last_status_time_ = this->now();
// }

// void BridgeNode::publish_ros_status_to_dds()
// {
//   std::optional<RosArray> copy;
//   {
//     std::lock_guard<std::mutex> lock(m_ros_);
//     if (!last_status_ros_.has_value()) return;
//     if (!status_dirty_) return;
//     copy = last_status_ros_;
//     status_dirty_ = false;
//   }

//   // 최소 크기 체크 (StatusType은 보통 18 기대)
//   if (copy->data.size() < 18) {
//     RCLCPP_WARN(get_logger(), "[ROS->DDS] StatusType ignored: data.size()=%zu", copy->data.size());
//     return;
//   }

//   const DdsStatus dds_msg = ros_to_dds_status(*copy);
//   dds_writer_status_.write(dds_msg);

//   RCLCPP_INFO(get_logger(),
//     "[ROS->DDS] StatusType write: USV(x,y,h)=(%.3f, %.3f, %.3f)",
//     copy->data[0], copy->data[1], copy->data[2]);
// }



void BridgeNode::poll_dds_and_publish_control()
{
  DdsControl latest;
  bool has = false;

  auto samples = dds_reader_control_.take();
  for (const auto& s : samples) {
    if (!s.info().valid()) continue;
    latest = s.data();
    has = true;
  }
  if (!has) return;

  RosArray ros = dds_to_ros_control(latest);
  ros_pub_control_->publish(ros);

  // ✅ echo loop 방지용 저장
  {
    std::lock_guard<std::mutex> lock(m_ros_);
    last_control_from_dds_ros_ = ros;
  }

  RCLCPP_INFO(get_logger(),
    "[DDS->ROS] USVControlType publish: rpm=%.0f steer=%.0f bow=%.0f",
    ros.data[0], ros.data[1], ros.data[2]);
}


void BridgeNode::poll_dds_and_publish_cmd_init()
{
  DdsCmdInit latest;
  bool has = false;

  // DDS에서 새 샘플들 가져오기 (take: 큐에서 꺼내기)
  auto samples = dds_reader_cmd_init_.take();
  for (const auto& s : samples)
  {
    if (!s.info().valid()) continue;
    latest = s.data();   // 가장 마지막(valid) 샘플을 최신으로 사용
    has = true;
  }
  if (!has) return;

  // DDS -> ROS 변환 후 publish
  RosArray ros = dds_to_ros_cmd_init(latest);
  ros_pub_cmd_init_->publish(ros);

  RCLCPP_INFO(get_logger(),
              "[DDS->ROS] CommandInitialType publish: commandID=%d",
              (int)latest.commandID());
}



} // namespace pils_dds_bridge


