from launch import LaunchDescription
from launch.actions import IncludeLaunchDescription
from launch.launch_description_sources import PythonLaunchDescriptionSource
from ament_index_python.packages import get_package_share_directory, get_package_prefix
from launch.actions import SetEnvironmentVariable


def generate_launch_description():
    return LaunchDescription([
        SetEnvironmentVariable(name='GAZEBO_PLUGIN_PATH',
                               value=[get_package_prefix('labamoro_plugin'), '/lib/labamoro_plugin']),
        SetEnvironmentVariable(name='GAZEBO_MODEL_PATH',
                               value=[get_package_share_directory('lab_amoro'), '/models']),
        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                [get_package_share_directory('gazebo_ros'), '/launch/gzserver.launch.py'])
        ),
        IncludeLaunchDescription(
            PythonLaunchDescriptionSource(
                [get_package_share_directory('gazebo_ros'), '/launch/gzclient.launch.py']),
            launch_arguments={'verbose': 'true'}.items()
        )
    ])
