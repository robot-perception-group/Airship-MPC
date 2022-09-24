#!/bin/bash

source ../../catkin_ws/devel/setup.bash

# stopping all rosnodes
rosnode kill --all
# stopping the gazebo client aka gui
killall gzclient
# stopping the gazebo server
killall gzserver
killall -9 roslaunch
killall -9 roslaunch
killall -9 roslaunch
killall server.py
killall ssd_server.sh
killall fw_simposix.elf
killall fc_server.sh

