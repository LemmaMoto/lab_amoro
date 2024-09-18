from lab_amoro.parallel_robot import *
from five_bar_models import *
from time import time
import sys


def main(args=None):
    # Initialize and start the ROS2 robot interface
    rclpy.init(args=args)
    robot = Robot("five_bar")
    start_robot(robot)

    # Get the error
    error_x = []
    error_y = []
    error_q12 = []
    error_q22 = []
    error_xD = []
    error_yD = []
    error_q12D = []
    error_q22D = []
    error_xDD = []
    error_yDD = []
    error_q12DD = []
    error_q22DD = []
    error_q11 = []
    error_q21 = []
    error_q11D = []
    error_q21D = []
    error_q11DD = []
    error_q21DD = []
    error_tau1 = []
    error_tau2 = []

    # Compare the models
    t_start = time()
    robot.start_oscillate()
    while time() < t_start + 60:
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
                # Error estimation
                error_y.append(abs(y - y_model))
                error_x.append(abs(x - x_model))
                error_q12.append(abs(q12 - q12_model))
                error_q22.append(abs(q22 - q22_model))

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
                # Error estimation
                error_yD.append(abs(yD - yD_model))
                error_xD.append(abs(xD - xD_model))
                error_q12D.append(abs(q12D - q12D_model))
                error_q22D.append(abs(q22D - q22D_model))

                # Test of the Direct Kinematic Model Second order
                # Data from gazebo
                q11DD = robot.active_left_joint.acceleration
                q21DD = robot.active_right_joint.acceleration
                xDD = robot.end_effector.acceleration_x
                yDD = robot.end_effector.acceleration_y
                q12DD = robot.passive_left_joint.acceleration
                q22DD = robot.passive_right_joint.acceleration
                # DKM 2nd order
                xDD_model, yDD_model = dkm2(q11, q12_model, q21, q22_model, q11D, q12D_model, q21D, q22D_model, q11DD,
                                            q21DD)
                q12DD_model, q22DD_model = dkm2_passive(q11, q12_model, q21, q22_model, q11D, q12D_model, q21D,
                                                        q22D_model, q11DD, q21DD, xDD_model, yDD_model)
                # Error estimation
                error_yDD.append(abs(yDD - yDD_model))
                error_xDD.append(abs(xDD - xDD_model))
                error_q12DD.append(abs(q12DD - q12DD_model))
                error_q22DD.append(abs(q22DD - q22DD_model))

                # Data from gazebo
                q11 = robot.active_left_joint.position
                q21 = robot.active_right_joint.position
                x = robot.end_effector.position_x
                y = robot.end_effector.position_y
                # Igm
                q11_model, q21_model = igm(x, y, -1, -1)
                # Error evaluation
                error_q11.append(abs(q11 - q11_model))
                error_q21.append(abs(q21 - q21_model))

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
                # Error evaluation
                error_q11D.append(abs(q11D - q11D_model))
                error_q21D.append(abs(q21D - q21D_model))

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
                # Error evaluation
                error_q11DD.append(abs(q11DD - q11DD_model))
                error_q21DD.append(abs(q21DD - q21DD_model))

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
                # Error evaluation
                error_tau1.append(abs(tau1 - tau_model[0]))
                error_tau2.append(abs(tau2 - tau_model[1]))

                # Go on with the oscillations
                robot.continue_oscillations()

        except KeyboardInterrupt:
            break

    # Print result and test pass
    print("Direct models:")
    for name, error in zip(
            ["x", "y", "q12", "q22", "xD", "yD", "q12D", "q22D", "xDD", "yDD", "q12DD", "q22DD"],
            [error_x, error_y, error_q12, error_q22, error_xD, error_yD, error_q12D, error_q22D, error_xDD, error_yDD,
             error_q12DD, error_q22DD]
    ):
        eval_text = "- OK" if sum(error) / len(error) < 1e-3 else "- error in model"
        print(name, " mean error: ", sum(error) / len(error), eval_text)
    print("Inverse models:")
    for name, error in zip(
            ["q11", "q21", "q11D", "q21D", "q11DD", "q21DD", "tau1", "tau2"],
            [error_q11, error_q21, error_q11D, error_q21D, error_q11DD, error_q21DD, error_tau1, error_tau2]
    ):
        eval_text = "- OK" if sum(error) / len(error) < 1e-3 else "- error in model"
        print(name, " mean error: ", sum(error) / len(error), eval_text)


if __name__ == "__main__":
    main(sys.argv)
