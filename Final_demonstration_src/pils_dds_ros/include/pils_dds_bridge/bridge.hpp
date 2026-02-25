#pragma once

#include <mutex>
#include <optional>
#include <string>

#include <rclcpp/rclcpp.hpp>
#include <std_msgs/msg/float64_multi_array.hpp>
#include <dds/dds.hpp>

#include "pils_dds_bridge/mapping.hpp"

namespace pils_dds_bridge
{
class BridgeNode : public rclcpp::Node
{
public:
  explicit BridgeNode(const rclcpp::NodeOptions& options = rclcpp::NodeOptions());

private:
  // params
  int domain_id_{200};
  double rate_hz_{10.0};

  // topic names (ROS == DDS)
  std::string topic_usv_control_{"USVControlType"};
  std::string topic_global_waypoint_{"GlobalWaypointType"};
  std::string topic_status_{"StatusType"};
  std::string topic_command_initial_{"CommandInitialType"};
  
  using RosArray = std_msgs::msg::Float64MultiArray;

  // ROS
  rclcpp::Subscription<RosArray>::SharedPtr ros_sub_control_; // ROS->DDS
  rclcpp::Publisher<RosArray>::SharedPtr ros_pub_gwp_;        // DDS->ROS
  rclcpp::Publisher<RosArray>::SharedPtr ros_pub_status_;     // DDS->ROS
  rclcpp::Publisher<RosArray>::SharedPtr ros_pub_cmd_init_;   // DDS->ROS

  std::mutex m_ros_;
  std::optional<RosArray> last_control_ros_;

  // DDS
  dds::domain::DomainParticipant participant_;
  dds::sub::Subscriber dds_sub_;
  dds::pub::Publisher dds_pub_;

  dds::topic::Topic<DdsControl> dds_topic_control_;
  dds::pub::DataWriter<DdsControl> dds_writer_control_;

  dds::topic::Topic<DdsGwp> dds_topic_gwp_;
  dds::sub::DataReader<DdsGwp> dds_reader_gwp_;

  dds::topic::Topic<DdsStatus> dds_topic_status_;
  dds::sub::DataReader<DdsStatus> dds_reader_status_;

  dds::topic::Topic<DdsCmdInit> dds_topic_cmd_init_;
  dds::sub::DataReader<DdsCmdInit> dds_reader_cmd_init_;


  rclcpp::TimerBase::SharedPtr timer_;

  bool control_dirty_{false};
  rclcpp::Time last_control_time_;


  rclcpp::Subscription<RosArray>::SharedPtr ros_sub_status_;   // ROS->DDS (StatusType)
  rclcpp::Publisher<RosArray>::SharedPtr ros_pub_control_;     // DDS->ROS (USVControlType)
  
  std::optional<RosArray> last_status_ros_;
  bool status_dirty_{false};
  rclcpp::Time last_status_time_;
  
  // Echo-loop guard (DDS->ROS로 방금 내보낸 메시지 저장)
  std::optional<RosArray> last_control_from_dds_ros_;
  std::optional<RosArray> last_status_from_dds_ros_;

  dds::sub::DataReader<DdsControl> dds_reader_control_;  // DDS->ROS (USVControlType)
  dds::pub::DataWriter<DdsStatus> dds_writer_status_;    // ROS->DDS (StatusType)
  void on_ros_status(const RosArray::SharedPtr msg);

  void publish_ros_status_to_dds();
  void poll_dds_and_publish_control();


  

  void init_parameters();
  void init_ros();
  void init_dds();

  void on_ros_control(const RosArray::SharedPtr msg);
  void tick_10hz();

  void publish_ros_control_to_dds();
  void poll_dds_and_publish_gwp();
  void poll_dds_and_publish_status();

  void poll_dds_and_publish_cmd_init();

};
} // namespace pils_dds_bridge
