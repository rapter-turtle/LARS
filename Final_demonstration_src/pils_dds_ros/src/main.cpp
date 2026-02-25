#include <rclcpp/rclcpp.hpp>
#include "pils_dds_bridge/bridge.hpp"

int main(int argc, char** argv)
{
  rclcpp::init(argc, argv);
  rclcpp::spin(std::make_shared<pils_dds_bridge::BridgeNode>());
  rclcpp::shutdown();
  return 0;
}
