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
    scope_x = Scope("Position x", 0.0, 0.2)
    scope_y = Scope("Position y", -0.5, 0.5)
    scope_q12 = Scope("Position q12", -1.0, 1.5)
    scope_q22 = Scope("Position q22", 0, 3.0)
    scope_xD = Scope("Velocity x", -0.02, 0.02)
    scope_yD = Scope("Velocity y", -0.02, 0.02)
    scope_q12D = Scope("Velocity q12", -0.2, 0.2)
    scope_q22D = Scope("Velocity q22", -0.2, 0.2)
    scope_xDD = Scope("Acceleration x", -0.05, 0.05)
    scope_yDD = Scope("Acceleration y", -0.05, 0.05)
    scope_q12DD = Scope("Acceleration q12", -0.5, 0.5)
    scope_q22DD = Scope("Acceleration q22", -0.5, 0.5)

    # Compare the models
    robot.start_oscillate()
    while True:  # Runs the simulation forever
        try:
            if robot.data_updated():
                # Test of the Direct Geometric Model
                # Data from gazebo
                q11 = robot.active_left_joint.position
                q21 = robot.active_right_joint.position
                x = robot.end_effector.position_x
                y = robot.end_effector.position_y
                q12 = robot.passive_left_joint.position
                q22 = robot.passive_right_joint.position
                # Dgm
                x_model, y_model = dgm(q11, q21, -1)
                q12_model, q22_model = dgm_passive(q11, q21, -1)
                # Plot update
                time = robot.get_time()
                scope_y.update(time, y, y_model)
                scope_x.update(time, x, x_model)
                scope_q12.update(time, q12, q12_model)
                scope_q22.update(time, q22, q22_model)

                # Test of the Direct Kinematic Model
                # Data from gazebo
                q11D = robot.active_left_joint.velocity
                q21D = robot.active_right_joint.velocity
                xD = robot.end_effector.velocity_x
                yD = robot.end_effector.velocity_y
                q12D = robot.passive_left_joint.velocity
                q22D = robot.passive_right_joint.velocity
                # Dkm
                xD_model, yD_model = dkm(q11, q12_model, q21, q22_model, q11D, q21D)
                q12D_model, q22D_model = dkm_passive(q11, q12_model, q21, q22_model, q11D, q21D, xD_model, yD_model)
                # Plot update
                time = robot.get_time()
                scope_yD.update(time, yD, yD_model)
                scope_xD.update(time, xD, xD_model)
                scope_q12D.update(time, q12D, q12D_model)
                scope_q22D.update(time, q22D, q22D_model)

                # Test of the Direct Kinematic Model Second order
                # Data from gazebo
                q11DD = robot.active_left_joint.acceleration
                q21DD = robot.active_right_joint.acceleration
                xDD = robot.end_effector.acceleration_x
                yDD = robot.end_effector.acceleration_y
                q12DD = robot.passive_left_joint.acceleration
                q22DD = robot.passive_right_joint.acceleration
                # DKM 2nd order
                xDD_model, yDD_model = dkm2(q11, q12_model, q21, q22_model, q11D, q12D_model, q21D, q22D_model, q11DD, q21DD)
                q12DD_model, q22DD_model = dkm2_passive(q11, q12_model, q21, q22_model, q11D, q12D_model, q21D, q22D_model, q11DD, q21DD, xDD_model, yDD_model)
                # Plot update
                time = robot.get_time()
                scope_yDD.update(time, yDD, yDD_model)
                scope_xDD.update(time, xDD, xDD_model)
                scope_q12DD.update(time, q12DD, q12DD_model)
                scope_q22DD.update(time, q22DD, q22DD_model)

                # Go on with the oscillations
                robot.continue_oscillations()

        except KeyboardInterrupt:
            break


if __name__ == "__main__":
    main(sys.argv)
