from lab_amoro.parallel_robot import *
from lab_amoro.plot_tools import *
from five_bar_models import *  # Modify this to test the biglide
import sys


def main(args=None):
    # Initialize and start the ROS2 robot interface
    rclpy.init(args=args)
    robot = Robot("five_bar")  # Modify this to test the biglide
    start_robot(robot)

    # Prepare plots - Comment scopes if you don't want to see them
    app = QtGui.QApplication([])
    scope_q11 = Scope("Position q11", -1.0, 1.0)
    scope_q21 = Scope("Position q21", -1.0, 1.0)
    scope_q11D = Scope("Velocity q11", -0.5, 0.5)
    scope_q21D = Scope("Velocity q21", -0.5, 0.5)
    scope_q11DD = Scope("Acceleration q11", -0.5, 0.5)
    scope_q21DD = Scope("Acceleration q21", -0.5, 0.5)
    scope_tau1 = Scope("Effort 1", -0.1, 0.1)
    scope_tau2 = Scope("Effort 2", -0.1, 0.1)

    # Compare the models
    robot.start_oscillate()
    try:
        while True:
            if robot.data_updated():
                # Test of the Inverse Geometric model
                # Data from gazebo
                q11 = robot.active_left_joint.position
                q21 = robot.active_right_joint.position
                x = robot.end_effector.position_x
                y = robot.end_effector.position_y
                # Igm
                q11_model, q21_model = igm(x, y, -1, -1)
                # Plot update
                time = robot.get_time()
                scope_q11.update(time, q11, q11_model)
                scope_q21.update(time, q21, q21_model)

                # Test of the Inverse Kinematic Model
                # Data from gazebo
                q12 = robot.passive_left_joint.position
                q22 = robot.passive_right_joint.position
                q11D = robot.active_left_joint.velocity
                q21D = robot.active_right_joint.velocity
                xD = robot.end_effector.velocity_x
                yD = robot.end_effector.velocity_y
                # Dkm
                q11D_model, q21D_model = ikm(q11, q12, q21, q22, xD, yD)
                # Plot update
                time = robot.get_time()
                scope_q11D.update(time, q11D, q11D_model)
                scope_q21D.update(time, q21D, q21D_model)

                # Test of the Inverse Kinematic Model Second order
                # Data from gazebo
                q12D = robot.passive_left_joint.velocity
                q22D = robot.passive_right_joint.velocity
                q11DD = robot.active_left_joint.acceleration
                q21DD = robot.active_right_joint.acceleration
                xDD = robot.end_effector.acceleration_x
                yDD = robot.end_effector.acceleration_y
                # IKM 2nd order
                q11DD_model, q21DD_model = ikm2(q11, q12, q21, q22, q11D, q12D, q21D, q22D, xDD, yDD)
                # Plot update
                time = robot.get_time()
                scope_q11DD.update(time, q11DD, q11DD_model)
                scope_q21DD.update(time, q21DD, q21DD_model)

                # Test of the Inverse Dynamic Model
                # Data from gazebo
                tau1 = robot.active_left_joint.effort
                tau2 = robot.active_right_joint.effort
                # Dynamic model
                M, c = dynamic_model(q11, q12, q21, q22, q11D, q12D, q21D, q22D)
                qDD = np.array([
                        q11DD,
                        q21DD
                    ])
                tau_model = M.dot(qDD) + c
                # Plot update
                scope_tau1.update(time, tau1, tau_model[0])
                scope_tau2.update(time, tau2, tau_model[1])

                robot.continue_oscillations()

    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main(sys.argv)
