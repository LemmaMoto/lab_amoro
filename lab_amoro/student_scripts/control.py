from lab_amoro.parallel_robot import *
from lab_amoro.plot_tools import *
from five_bar_models import *  # Modify this to use the biglide
import sys


def main(args=None):
    # Initialize and start the ROS2 robot interface
    rclpy.init(args=args)
    robot = Robot("five_bar")  # Modify this to use the biglide
    start_robot(robot)

    # Prepare plots
    app = QtGui.QApplication([])
    scope_joint1 = Scope("Joint 1", -0.5, 1.5)
    scope_joint2 = Scope("Joint 2", -1.5, 1.5)

    # Create the trajectory as arrays in Cartesian space (position, velocity, acceleration)

    # Create the trajectory as arrays in joint space using the inverse models (position, velocity, acceleration)

    index = 0

    # Controller
    try:
        robot.apply_efforts(0.0, 0.0)  # Required to start the simulation
        while True:
            if robot.data_updated():
                # Robot available data - This is the only data thet you can get from a real robot (joint encoders)
                q11 = robot.active_left_joint.position
                q21 = robot.active_right_joint.position
                q11D = robot.active_left_joint.velocity
                q21D = robot.active_right_joint.velocity

                # CTC controller
                # TODO
                robot.apply_efforts(tau_left, tau_right)

                # Scope update
                time = robot.get_time()
                if time < 5.0:
                    scope_joint1.update(time, 0.0, 0.0)
                    scope_joint2.update(time, 0.0., 0.0)

                if index < len(trajectory)-1:
                    index += 1  # Next point in trajectory

    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main(sys.argv)
