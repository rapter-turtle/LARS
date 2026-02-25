cmake_minimum_required(VERSION 3.8)
project(pils_dds_ros)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_options(-Wall -Wextra -Wpedantic)

# =========================
# ROS2
# =========================
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)
find_package(GeographicLib REQUIRED)

# =========================
# RTI Connext DDS
# =========================
if(NOT DEFINED ENV{NDDSHOME})
  message(FATAL_ERROR "NDDSHOME environment variable is not set. Please set it to RTI Connext DDS installation directory.")
endif()
set(NDDSHOME $ENV{NDDSHOME})

# 너 환경에 맞게 colcon build --cmake-args -DDDS_TARGET_ARCH=... 로 바꿀 수 있게 캐시 처리
set(DDS_TARGET_ARCH "x64Linux4gcc7.3.0" CACHE STRING "RTI DDS target arch folder under NDDSHOME/lib")

# =========================
# Include / Link dirs
# =========================
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/dds_gen

  ${NDDSHOME}/include
  ${NDDSHOME}/include/ndds
  ${NDDSHOME}/include/ndds/hpp
)

link_directories(
  ${NDDSHOME}/lib/${DDS_TARGET_ARCH}
)

# =========================
# DDS compile defs
# =========================
set(DDS_COMPILE_DEFS
  RTI_UNIX
  RTI_LINUX
  RTI_64BIT
  RTI_STATIC
)

# =========================
# Bridge node
# =========================
add_executable(pils_dds_bridge_node
  src/main.cpp
  src/bridge.cpp
  dds_gen/pils.cxx
  dds_gen/pilsPlugin.cxx
  # 만약 rtiddsgen이 Support 파일도 만들면 필요할 수 있음:
  # dds_gen/pilsSupport.cxx
)

target_compile_definitions(pils_dds_bridge_node PRIVATE ${DDS_COMPILE_DEFS})

ament_target_dependencies(pils_dds_bridge_node
  rclcpp
  std_msgs
)

target_link_libraries(pils_dds_bridge_node
  # RTI DDS static libs (너가 기존에 쓰던 z 붙은 버전 유지)
  nddscpp2z
  nddscz
  nddscorez
  # system
  m
  pthread
  dl
  rt
  GeographicLib::GeographicLib
)

install(TARGETS pils_dds_bridge_node
  DESTINATION lib/${PROJECT_NAME}
)


ament_package()
