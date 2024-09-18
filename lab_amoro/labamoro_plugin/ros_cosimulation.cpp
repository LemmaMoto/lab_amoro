//
// Created by Damien SIX on 06/08/2021.
//

#include "ros_cosimulation.h"
#include <gazebo/physics/Model.hh>
#include <gazebo/physics/Joint.hh>
#include <gazebo/physics/Link.hh>
#include <gazebo/physics/World.hh>

namespace gazebo
{
    GZ_REGISTER_MODEL_PLUGIN(GazeboRosCosimulation)

    GazeboRosCosimulation::GazeboRosCosimulation()
    : ModelPlugin(), effortLeft(0.0), effortRight(0.0),logger(rclcpp::get_logger("gazebo_cosimulation"))
    {}

    GazeboRosCosimulation::~GazeboRosCosimulation() = default;

    void GazeboRosCosimulation::Load(gazebo::physics::ModelPtr model, sdf::ElementPtr sdf) {
        RCLCPP_INFO(logger, "Starting cosimulation plugin");

        world = model->GetWorld();
        world->SetPaused(true);
        nextStopTime = 0.0;

        active_left_joint = model->GetJoint("A11");
        active_right_joint = model->GetJoint("A21");
        passive_left_joint = model->GetJoint("A12");
        passive_right_joint = model->GetJoint("A22");
        end_effector = model->GetLink("endEffector");
        base = model->GetLink("base");

        ros_node_ = gazebo_ros::Node::Get();
        joint_efforts_sub = ros_node_->create_subscription<lab_amoro_interfaces::msg::JointEfforts>(
                "joint_efforts",
                1,
                std::bind(&GazeboRosCosimulation::OnRosJointEfforts, this, std::placeholders::_1)
                );
        run_sub = ros_node_->create_subscription<std_msgs::msg::Float32>(
                "simulation_run_until",
                1,
                std::bind(&GazeboRosCosimulation::RunUntil, this, std::placeholders::_1)
        );
        robot_states_pub = ros_node_->create_publisher<lab_amoro_interfaces::msg::RobotStates>(
                "robot_states",
                1
        );
        updateConnection = event::Events::ConnectWorldUpdateBegin(
                std::bind(&GazeboRosCosimulation::OnUpdate, this));
        readyForNewStep = true;
    }

    void GazeboRosCosimulation::OnUpdate() {
        active_left_joint->SetForce(0, effortLeft);
        active_right_joint->SetForce(0, effortRight);
        if (world->SimTime().Float() >= nextStopTime && !world->IsPaused()) {
            world->SetPaused(true);
            // Publish states
            lab_amoro_interfaces::msg::RobotStates states;
            states.end_effector_x.position = end_effector->WorldPose().X() - base->WorldPose().X();
            states.end_effector_x.velocity = end_effector->WorldLinearVel().X();
            states.end_effector_x.acceleration = end_effector->WorldLinearAccel().X();
            states.end_effector_y.position = end_effector->WorldPose().Y() - base->WorldPose().Y();
            states.end_effector_y.velocity = end_effector->WorldLinearVel().Y();
            states.end_effector_y.acceleration = end_effector->WorldLinearAccel().Y();
            states.active_left_joint.position = active_left_joint->Position();
            states.active_left_joint.velocity = active_left_joint->GetVelocity(0);
            states.active_right_joint.position = active_right_joint->Position();
            states.active_right_joint.velocity = active_right_joint->GetVelocity(0);
            states.passive_left_joint.position = passive_left_joint->Position();
            states.passive_left_joint.velocity = passive_left_joint->GetVelocity(0);
            states.passive_right_joint.position = passive_right_joint->Position();
            states.passive_right_joint.velocity = passive_right_joint->GetVelocity(0);
            states.simulation_time = world->SimTime().Float();
            robot_states_pub->publish(states);
            readyForNewStep = true;
        }
    }

    void GazeboRosCosimulation::OnRosJointEfforts(const lab_amoro_interfaces::msg::JointEfforts::SharedPtr msg) {
        effortLeft = msg->left_joint;
        effortRight = msg->right_joint;
    }

    void GazeboRosCosimulation::RunUntil(const std_msgs::msg::Float32::SharedPtr msg){
        // We accept a new step only if we are not already running one
        if (readyForNewStep) {
            nextStopTime = msg->data;
            readyForNewStep = false;
            world->SetPaused(false);
        }
    }

}