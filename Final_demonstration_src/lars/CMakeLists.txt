cmake_minimum_required(VERSION 3.8)
project(lars)

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(aura_msg REQUIRED)   # ⬅ 추가!
 
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ACADOS PATH
set(ACADOS_ROOT "/home/user/acados")

include_directories(
    ${ACADOS_ROOT}/include
    ${ACADOS_ROOT}/include/acados
    ${ACADOS_ROOT}/include/blasfeo/include
    ${ACADOS_ROOT}/include/hpipm/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

link_directories(
    ${ACADOS_ROOT}/lib
)

# ===============================================================
# Gather all acados generated C code
# ===============================================================
# file(GLOB GEN_C_SOURCES
#     "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/*.c"
#     "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/*/*.c"
# )
file(GLOB GEN_C_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/acados_solver_heron.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/acados_sim_solver_heron.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_model/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_constraints/*.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_cost/*.c"
)





# Remove auto-generated mains (acados generates its own test main)
list(FILTER GEN_C_SOURCES EXCLUDE REGEX ".*main_.*\\.c$")

# ===============================================================
# Build executable
# ===============================================================
# add_executable(heron_mpc
#     src/main_heron_mpc.cpp
#     ${GEN_C_SOURCES}
# )
add_executable(heron_mpc_module
    src/main_heron_module.cpp
    ${GEN_C_SOURCES}
)

add_executable(heron_los
    src/main_heron_los.cpp
    ${GEN_C_SOURCES}
)

# add_executable(heron_mpc
#     src/main_heron_mpc.cpp
#     ${GEN_C_SOURCES}
# )
# # target_include_directories(heron_mpc PRIVATE
# #     ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code
# # )
# target_include_directories(heron_mpc PRIVATE
#     ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code
#     ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_model
#     ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_constraints
#     ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_cost
# )



# target_link_libraries(heron_mpc
#     acados
#     hpipm
#     blasfeo
#     m
#     pthread
# )

# ament_target_dependencies(heron_mpc
#     rclcpp
#     aura_msg      # ⬅ 추가!
# )


# install(TARGETS heron_mpc
#     DESTINATION lib/${PROJECT_NAME}
# )


target_include_directories(heron_mpc_module PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_model
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_constraints
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_cost
)



target_link_libraries(heron_mpc_module
    acados
    hpipm
    blasfeo
    m
    pthread
)

ament_target_dependencies(heron_mpc_module
    rclcpp
    aura_msg      # ⬅ 추가!
)


install(TARGETS heron_mpc_module
    DESTINATION lib/${PROJECT_NAME}
)


target_include_directories(heron_los PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_model
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_constraints
    ${CMAKE_CURRENT_SOURCE_DIR}/c_generated_code/heron_cost
)



target_link_libraries(heron_los
    acados
    hpipm
    blasfeo
    m
    pthread
)

ament_target_dependencies(heron_los
    rclcpp
    aura_msg      # ⬅ 추가!
)


install(TARGETS heron_los
    DESTINATION lib/${PROJECT_NAME}
)


install(
    DIRECTORY c_generated_code/
    DESTINATION share/${PROJECT_NAME}/c_generated_code
)

ament_package()
