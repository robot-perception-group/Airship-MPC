# Airship-MPC
![airship_cover_with_inlay](https://user-images.githubusercontent.com/32105268/192122387-465e4489-5635-44d2-ab37-2071c48651a0.png)

Perception-driven Formation Control of Airships for project AirCap using nonlinear model predictive control (MPC)

A preprint describing our method is available on https://arxiv.org/abs/2209.13040

image:: https://img.youtube.com/vi/ihS0_VRD_kk/maxresdefault.jpg
    :alt: Video: Perception-driven Formation Control of Airships
    :target: https://www.youtube.com/watch?v=ihS0_VRD_kk


## Usage

1. Prepare a system matching the pre-requirements
2. Download this sourcecode repository into your home folder
3. Install the compilation requirements and compile code
4. Run the simulation demo or the MPC evaluation demos

### Pre Requirements

This code has been tested with Ubuntu Linux 20.04 LTS (Focal Fossa) and ROS Noetic Ninjemys only.
The Neural Network detector used for person detection requires an NVIDIA GPU with current Cuda capable drivers installed.

### Download

You can download the code with the command  ```git clone https://github.com/robot-perception-group/Airship-MPC```

### Compile

Please run the script ```./install_and_compile.sh``` it will install all necessary requirements including ROS Noetic Ninjemys and compile the code.

### Execute

 * The Gazebo simulation with wind and person tracking can be started with ```cd experiments/sim; ./airship_sim.sh 3```
 * The MPC can be tested standalone using it's own motion model for simulation with ```cd src; python3 blimps_costtest.py```
 * The solver generated trajectories and plots from our paper can be recreated with the scripts found in the subfolders in folder experiments/trajectory

### Modify
 * A simulation with a stationary person can be started with <br> ```cd experiments/sim; ./airship_sim.sh 3 100 test arena_BLIMP_stat_target```
 * Constraints and properties of the MPC can be modified in ```src/blimp_nmpc_wrapper_node/nodes/formation_config.py``` including the number of airships in the formation. Whenever this is changed, a recompilation is necessary.
 * ```submodules/AirCap/packages/3rdparty/airship_simulation/deflate_blimp.sh``` can be used to alter the airship buoyancy and rigidity
 * Wind and other environmental parameters affecting the airships can be modified by editing ```submodules/AirCap/packages/3rdparty/airship_simulation/blimp_description/urdf/description_plugin.xacro```. This folder also includes the physical model description of the airship and it's aerodynamic properties in Gazebo URDF format.
 * The scripts in ```submodules/AirCap/scripts/simulation/``` can be used to run various other experiments, including standard AirCap with multicopters

### Real Airship
To run the controller on a real airship, the vehicle needs an OpenPilot Revolution Flight controller or compatible
and a companion computer capable of running neural networks and with wireless communication capability connected via USB. We are using an NVIDIA Jetson TX1 or TX2 with integrated wifi.
ROS should ideally be run in a multi master setup using fkie_multimaster, with an instance of the MPC running on each vehicle.
