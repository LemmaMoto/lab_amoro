//
// Created by Damien SIX on 06/08/2021.
//

#ifndef LAB_AMORO_ROS_COSIMULATION_H
#define LAB_AMORO_ROS_COSIMULATION_H
#include <gazebo/common/Plugin.hh>
#include <lab_amoro_interfaces/msg/joint_efforts.hpp>
#include <lab_amoro_interfaces/msg/robot_states.hpp>
#include <std_msgs/msg/float32.hpp>
#include <rclcpp/rclcpp.hpp>
#include <gazebo_ros/node.hpp>

namespace gazebo
{

    class GazeboRosCosimulation : public gazebo::ModelPlugin
            {
            public:
                /// Constructor
                GazeboRosCosimulation();

                /// Destructor
                virtual ~GazeboRosCosimulation();

            protected:
                void Load(gazebo::physics::ModelPtr model, sdf::ElementPtr sdf) override;

                /// Optional callback to be called at every simulation iteration.
                virtual void OnUpdate();

            private:
                /// ROS2 Callbacks
                void OnRosJointEfforts(lab_amoro_interfaces::msg::JointEfforts::SharedPtr msg);
                void RunUntil(std_msgs::msg::Float32::SharedPtr msg);

                /// ROS2 Node
                gazebo_ros::Node::SharedPtr ros_node_;

                float effortLeft;
                float effortRight;
                float nextStopTime;
                bool readyForNewStep;
                gazebo::physics::WorldPtr world;
                gazebo::physics::JointPtr active_left_joint;
                gazebo::physics::JointPtr active_right_joint;
                gazebo::physics::JointPtr passive_left_joint;
                gazebo::physics::JointPtr passive_right_joint;
                gazebo::physics::LinkPtr end_effector;
                gazebo::physics::LinkPtr base;
                rclcpp::Subscription<lab_amoro_interfaces::msg::JointEfforts>::SharedPtr joint_efforts_sub;
                rclcpp::Subscription<std_msgs::msg::Float32>::SharedPtr run_sub;
                rclcpp::Publisher<lab_amoro_interfaces::msg::RobotStates>::SharedPtr robot_states_pub;

                rclcpp::Logger logger;

                private: event::ConnectionPtr updateConnection;
            };

}
#endif //LAB_AMORO_ROS_COSIMULATION_H
