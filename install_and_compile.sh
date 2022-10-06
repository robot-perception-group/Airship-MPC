#!/bin/bash

echo "This will install all dependencies for Airship MPC simulation and start the compilation."
echo "Press Enter to continue, CTRL+C to abort."
echo "Warning, this will call sudo to install packages with apt-get."
echo "If in doubt, read this file and execute the commands one by one with changes as necessary."
read disclaimer

echo "Will install required system packages with apt-get"
sudo apt-get update
sudo apt-get install git vim lsb-release build-essential curl llvm-dev libclang-dev clang python-is-python3 python3-pip libgoogle-glog-dev python2 screen psmisc bc || { echo "Failure"; exit;}


echo "Checking for ROS Noetic"
if [ ! -e /opt/ros/noetic ]; then
	echo "No ROS Noetic found, will install..."
	sudo sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list'
	curl -s https://raw.githubusercontent.com/ros/rosdistro/master/ros.asc | sudo apt-key add -
	sudo apt-get update;
else
	echo "ROS Noetic is already installed in /opt, proceeding";
fi

echo "Will install required ROS packages"
sudo apt-get install ros-noetic-ros-base gazebo11 ros-noetic-gazebo-dev ros-noetic-octomap-msgs ros-noetic-tf ros-noetic-mrpt2 ros-noetic-image-transport ros-noetic-rviz ros-noetic-cv-camera ros-noetic-image-geometry ros-noetic-gazebo-plugins ros-noetic-octomap-ros ros-noetic-controller-manager ros-noetic-joint-state-publisher ros-noetic-robot-state-publisher ros-noetic-xacro ros-noetic-gazebo-ros-control ros-noetic-position-controllers ros-noetic-compressed-image-transport || { echo "Failure"; exit;}
source /opt/ros/noetic/setup.bash

echo "Will install required GIT submodules"
git submodule update --init --recursive || { echo "Failure"; exit;}

echo "Will install python modules"
python3 -m pip install -U pip
python3 -m pip install opengen torch torchvision torchaudio matplotlib
export PATH="$HOME/.local/bin:$PATH"

echo "Will install rust"
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh /dev/stdin -y
export PATH="$HOME/.cargo/bin:$PATH"

echo "Will compile MPC solver into ROS module"
cd src
python3 blimps_ros.py || { echo "Failure"; exit;}
cd ..

echo "Will compile ROS packages"
mkdir catkin_ws
mkdir catkin_ws/src
cd catkin_ws/src
catkin_init_workspace
ln -s ../../submodules
ln -s ../../src
cd ..
touch src/submodules/AirCap/packages/optional/basler_image_capture/CATKIN_IGNORE
touch src/submodules/AirCap/packages/optional/ptgrey_image_capture/CATKIN_IGNORE
catkin_make || { echo "Failure"; exit;}
source devel/setup.bash
cd ..

echo "Will compile Librepilot Firmware for SITL simulation"
cd submodules/AirCap/packages/3rdparty/airship_simulation/LibrePilot/
make simposix || { echo "Failure"; exit;}
cd ../../../../../..

echo -e "\n\n Success!"
echo "You can run the GAZEBO simulation with the following commands:"
echo "  cd experiments/sim"
echo "  ./airship_sim.sh 3"
echo "You can ithen stop the simulation with the command"
echo "  ./cleanup.sh"


