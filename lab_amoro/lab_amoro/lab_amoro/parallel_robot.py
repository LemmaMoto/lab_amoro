import rclpy
from threading import Thread
from rclpy.node import Node
from std_msgs.msg import Float32
from std_srvs.srv import Empty
from lab_amoro_interfaces.msg import *


class Joint:
    def __init__(self, joint_name, joint_offset=0.0):
        self.joint_name = joint_name
        self.offset = joint_offset

    position = 0.0
    velocity = 0.0
    acceleration = 0.0
    effort = 0.0


class EndEffector:
    position_x = 0.0
    position_y = 0.0
    velocity_x = 0.0
    velocity_y = 0.0
    acceleration_x = 0.0
    acceleration_y = 0.0


class Robot(Node):
    def __init__(self, robot_type):
        super().__init__('robot_interface')
        # Oscillations
        self.oscillate = False
        if robot_type == 'five_bar':
            self.oscillation_effort = 0.0002
            self.oscillation_limit = 0.05
        else:
            self.oscillation_limit = 0.01
            self.oscillation_effort = 0.01
        self.oscillate_effort_left = self.oscillation_effort
        self.oscillate_effort_right = -self.oscillation_effort
        # Simulation step topic
        self.simulation_step_pub = self.create_publisher(Float32,
                                                          'simulation_run_until',
                                                          1)
        # Joint state topics
        self.create_subscription(RobotStates,
                                 'robot_states',
                                 self.robot_states_callback,
                                 1)
        # Joint effort topic
        self.joint_efforts_pub = self.create_publisher(JointEfforts,
                                                       'joint_efforts',
                                                       1)
        # Joints
        self.active_left_joint = Joint("A11")
        self.active_right_joint = Joint("A21")
        if robot_type == "five_bar":
            self.passive_left_joint = Joint("A12", 0.8558782490506787)
            self.passive_right_joint = Joint("A22", 2.2857144045391147)
        else:
            self.passive_left_joint = Joint("A12", 0.7853981192930068)
            self.passive_right_joint = Joint("A22", 2.3561945342967863)
        # EndEffector
        self.end_effector = EndEffector()
        # Reset simulation
        self.reset_client = self.create_client(Empty, 'reset_simulation')
        self.reset_simulation()
        # Update tracking
        self.updated = False
        # Previous time for acceleration in joint computation
        self.previous_sim_time = 0.0
        # Time for next step
        self.next_step = 0.0

    def robot_states_callback(self, msg):
        dt = msg.simulation_time - self.previous_sim_time
        if not dt == 0.0:
            self.previous_sim_time = msg.simulation_time
            self.active_left_joint.acceleration = (msg.active_left_joint.velocity - self.active_left_joint.velocity) / dt
            self.passive_left_joint.acceleration = (msg.passive_left_joint.velocity - self.passive_left_joint.velocity) / dt
            self.active_right_joint.acceleration = (msg.active_right_joint.velocity - self.active_right_joint.velocity) / dt
            self.passive_right_joint.acceleration = (msg.passive_right_joint.velocity - self.passive_right_joint.velocity) / dt
            self.active_left_joint.position = msg.active_left_joint.position + self.active_left_joint.offset
            self.passive_left_joint.position = msg.passive_left_joint.position + self.passive_left_joint.offset
            self.active_right_joint.position = msg.active_right_joint.position + self.active_right_joint.offset
            self.passive_right_joint.position = msg.passive_right_joint.position + self.passive_right_joint.offset
            self.active_left_joint.velocity = msg.active_left_joint.velocity
            self.passive_left_joint.velocity = msg.passive_left_joint.velocity
            self.active_right_joint.velocity = msg.active_right_joint.velocity
            self.passive_right_joint.velocity = msg.passive_right_joint.velocity
            self.end_effector.acceleration_x = msg.end_effector_x.acceleration
            self.end_effector.acceleration_y = msg.end_effector_y.acceleration
            self.end_effector.position_x = msg.end_effector_x.position
            self.end_effector.position_y = msg.end_effector_y.position
            self.end_effector.velocity_x = msg.end_effector_x.velocity
            self.end_effector.velocity_y = msg.end_effector_y.velocity
        self.updated = True

    def continue_oscillations(self):
        if self.oscillate:
            if self.oscillate_effort_left > 0.0 and self.active_left_joint.position > self.oscillation_limit:
                self.oscillate_effort_left = -self.oscillation_effort
            if self.oscillate_effort_left < 0.0 and self.active_left_joint.position < -self.oscillation_limit:
                self.oscillate_effort_left = self.oscillation_effort
            if self.oscillate_effort_right > 0.0 and self.active_right_joint.position > self.oscillation_limit:
                self.oscillate_effort_right = -self.oscillation_effort
            if self.oscillate_effort_right < 0.0 and self.active_right_joint.position < -self.oscillation_limit:
                self.oscillate_effort_right = self.oscillation_effort
            self.apply_efforts(self.oscillate_effort_left, self.oscillate_effort_right)

    def reset_simulation(self):
        service_available = self.reset_client.wait_for_service(10)
        if service_available:
            self.reset_client.call_async(Empty.Request())
            self.next_step = 0.0
        else:
            self.get_logger().warning("Gazebo reset service not available. Please start Gazebo first")

    def start_oscillate(self):
        self.oscillate = True
        self.apply_efforts(0.0, 0.0)

    def stop_oscillate(self):
        self.oscillate = False
        self.reset_simulation()

    def data_updated(self):
        updated = self.updated
        if updated:
            self.updated = False
        return updated

    def get_time(self):
        return self.previous_sim_time

    def apply_efforts(self, left_joint, right_joint):
        efforts = JointEfforts()
        efforts.right_joint = right_joint
        efforts.left_joint = left_joint
        self.active_left_joint.effort = left_joint
        self.active_right_joint.effort = right_joint
        self.joint_efforts_pub.publish(efforts)
        next_step = Float32()
        self.next_step += 0.01  # 10 ms timestep
        next_step.data = self.next_step
        self.simulation_step_pub.publish(next_step)


def start_robot(node):
    node_thread = Thread(target=rclpy.spin, args=(node,))
    node_thread.daemon = True
    node_thread.start()
