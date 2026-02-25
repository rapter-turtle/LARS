#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <std_msgs/msg/float64_multi_array.hpp>
#include "pils_dds_bridge/uvr_filter.hpp" 
// rtiddsgen 생성 헤더 (이름 다르면 수정)
#include "pils.hpp"
#include <GeographicLib/UTMUPS.hpp>
#include <chrono>  

struct UTMResult {
  double easting;
  double northing;
  int zone;
  bool northp;
};

inline UTMResult ll_to_utm(double lat_deg, double lon_deg)
{
  UTMResult r{};
  GeographicLib::UTMUPS::Forward(lat_deg, lon_deg, r.zone, r.northp, r.easting, r.northing);
  return r;
}

namespace pils_dds_bridge
{
using RosArray = std_msgs::msg::Float64MultiArray;

// DDS 타입 alias
using DdsControl = ULARS::EO::USVControlType;
using DdsGwp     = ULARS::C2::GlobalWaypointType;
using DdsWp      = ULARS::C2::Waypoint;
using DdsStatus  = ULARS::C2::StatusType;
using DdsCmdInit = ULARS::C2::CommandInitialType;

inline short clamp_to_short(double v)
{
  if (std::isnan(v) || std::isinf(v)) return 0;
  v = std::round(v);
  v = std::max<double>(v, static_cast<double>(std::numeric_limits<short>::min()));
  v = std::min<double>(v, static_cast<double>(std::numeric_limits<short>::max()));
  return static_cast<short>(v);
}

// ---------------------------
// USVControlType: ROS -> DDS
// array: [RPM, Steering, BowThrust]
// ---------------------------
inline DdsControl ros_to_dds_control(const RosArray& ros)
{
  DdsControl dds;
  const auto& a = ros.data;

  const double rpm   = (a.size() > 0) ? a[0] : 0.0;
  const double steer = (a.size() > 1) ? a[1] : 0.0;
  const double bow   = (a.size() > 2) ? a[2] : 0.0;

  dds.integratedTargetRPM(clamp_to_short(rpm));
  dds.integratedTargetSteering(clamp_to_short(steer));
  dds.integratedBowThrust(clamp_to_short(bow));
  return dds;
}

// ---------------------------
// GlobalWaypointType: DDS -> ROS
// array: [commandID, N, (lat,lon,lateralTol,goalTol,maxVel)*N]
// ---------------------------
inline RosArray dds_to_ros_gwp(const DdsGwp& dds)
{
  RosArray ros;
  const double cmd = static_cast<double>(dds.commandID());
  const auto& seq = dds.waypointData();
  const std::size_t N = seq.size();

  ros.data.reserve(2 + 5 * N);
  ros.data.push_back(cmd);
  ros.data.push_back(static_cast<double>(N));

  for (std::size_t i = 0; i < N; ++i)
  {
    const DdsWp& wp = seq[i];
    ros.data.push_back(wp.latitude());
    ros.data.push_back(wp.longitude());
    ros.data.push_back(static_cast<double>(wp.lateralErrorTolerance()));
    ros.data.push_back(static_cast<double>(wp.goalErrorTolerance()));
    ros.data.push_back(static_cast<double>(wp.maxVelocity()));
  }
  return ros;
}

// ---------------------------
// StatusType: DDS -> ROS (15개 IDL 순서)
// ---------------------------
inline RosArray dds_to_ros_status(const DdsStatus& s)
{
  RosArray ros;
  ros.data.resize(18);

  // const double lat = s.USV_x();
  // const double lon = s.USV_y();
  // const auto utm = ll_to_utm(lat, lon);
  // static double t_sec = 0.0;
  // t_sec += 0.1;

  // static uvr_filter::PoseLPFEKF pose_est_;
  // uvr_filter::PoseLPFEKF::Output out;
  // if (!pose_est_.updatePose(utm.easting, utm.northing, uvr_filter::wrapToPi(-s.USV_h() + 0.5*3.141592), t_sec, out)) {
  //     // Handle error if necessary
  //     return ros;
  // }

  double usv_x = s.USV_x();
  double craddle_x1 = s.Craddle_x1();
  double craddle_x2 = s.Craddle_x2();
  double craddle_x3 = s.Craddle_x3();
  double craddle_x4 = s.Craddle_x4();
  double mship_x = s.MSHIP_x();

  double usv_h_offset = 0.0;




  // -------- 실제 시간 기반 t_sec 누적 --------
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

  // dt가 비정상적으로 크거나 작으면(일시정지/디버깅 등) 클램프 권장
  if (dt < 0.0) dt = 0.0;
  if (dt > 0.5) dt = 0.5;   // 예: 최대 0.5s로 제한(원하면 조정)

  t_sec += dt;
  // ------------------------------------------


  // static double t_sec = 0.0;
  // t_sec += 0.1;

  static uvr_filter::PoseLPFEKF pose_est_;
  uvr_filter::PoseLPFEKF::Output out;
  // if (!pose_est_.updatePose(usv_x, s.USV_y(), uvr_filter::wrapToPi(-s.USV_h() + usv_h_offset), t_sec, out)) {
  //     // Handle error if necessary
  //     return ros;
  // }


  if (!pose_est_.updatePose(usv_x, s.USV_y(), uvr_filter::wrapToPi(s.USV_h() + usv_h_offset), t_sec, out)) {
    // Handle error if necessary
    return ros;
  }
  // std::cout << "u : " << out.u << ", v : " << out.v << ", r : " << out.r << std::endl;


  ros.data[0]  = usv_x;
  ros.data[1]  = s.USV_y();
  ros.data[2]  = s.USV_h();
  ros.data[3]  = out.u;
  ros.data[4]  = out.v;
  ros.data[5]  = out.r;
  ros.data[6]  = mship_x;
  ros.data[7]  = s.MSHIP_y();
  ros.data[8]  = s.MSHIP_h();
  ros.data[9]  = craddle_x1;
  ros.data[10]  = s.Craddle_y1();
  ros.data[11]  = craddle_x2;
  ros.data[12]  = s.Craddle_y2();
  ros.data[13] = craddle_x3;
  ros.data[14] = s.Craddle_y3();
  ros.data[15] = craddle_x4;
  ros.data[16] = s.Craddle_y4();
  ros.data[17] = s.Craddle_h();

  return ros;
}



// ---------------------------
// USVControlType: DDS -> ROS
// array: [RPM, Steering, BowThrust]
// ---------------------------
inline RosArray dds_to_ros_control(const DdsControl& dds)
{
  RosArray ros;
  ros.data.resize(3);
  ros.data[0] = static_cast<double>(dds.integratedTargetRPM());
  ros.data[1] = static_cast<double>(dds.integratedTargetSteering());
  ros.data[2] = static_cast<double>(dds.integratedBowThrust());
  return ros;
}

// ---------------------------
// StatusType: ROS -> DDS
// ROS array (same as dds_to_ros_status output):
// [0]usv_x [1]usv_y [2]usv_h [3]u [4]v [5]r
// [6]mship_x [7]mship_y [8]mship_h
// [9]cx1 [10]cy1 [11]cx2 [12]cy2 [13]cx3 [14]cy3 [15]cx4 [16]cy4
// [17]craddle_h
// ---------------------------
inline DdsStatus ros_to_dds_status(const RosArray& ros)
{
  DdsStatus s;
  const auto& a = ros.data;
  auto get = [&](std::size_t i) -> double { return (i < a.size()) ? a[i] : 0.0; };

  s.USV_x(get(0));
  s.USV_y(get(1));
  s.USV_h(get(2));

  // Inverse of your DDS->ROS x transforms
  s.MSHIP_x(get(6));
  s.MSHIP_y(get(7));
  s.MSHIP_h(get(8));

  s.Craddle_x1(get(9));
  s.Craddle_y1(get(10));
  s.Craddle_x2(get(11));
  s.Craddle_y2(get(12));
  s.Craddle_x3(get(13));
  s.Craddle_y3(get(14));
  s.Craddle_x4(get(15));
  s.Craddle_y4(get(16));
  s.Craddle_h(get(17));

  return s;
}

// ---------------------------
// CommandInitialType: DDS -> ROS
// ROS array: [commandID]
// ---------------------------
inline RosArray dds_to_ros_cmd_init(const DdsCmdInit& dds)
{
  RosArray ros;
  ros.data.resize(1);
  ros.data[0] = static_cast<double>(dds.commandID());
  return ros;
}


} // namespace pils_dds_bridge


