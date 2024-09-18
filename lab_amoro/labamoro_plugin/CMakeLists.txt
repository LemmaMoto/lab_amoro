cmake_minimum_required(VERSION 3.5)
project(labamoro_plugin)

find_package(ament_cmake REQUIRED)
find_package(lab_amoro_interfaces REQUIRED)
find_package(std_msgs REQUIRED)
find_package(rclcpp REQUIRED)
find_package(gazebo_ros REQUIRED)

include_directories(include)

add_library(gazebo_ros_cosimulation SHARED ros_cosimulation.cpp)
ament_target_dependencies(gazebo_ros_cosimulation "lab_amoro_interfaces" "std_msgs" "rclcpp" "gazebo_ros")

install(TARGETS gazebo_ros_cosimulation DESTINATION lib/${PROJECT_NAME})

ament_package()