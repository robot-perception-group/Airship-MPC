# Airship-MPC

Nonlinear Model Predictive Formation Controller for Airships for project AirCap.

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
 * The solver generated trajectories and plots from our paper can be recreated with the scripts found in the folder experiments/trajectory

### Modify
 * A simulation with a stationary person can be started with <br> ```cd experiments/sim; ./airship_sim.sh 3 100 test arena_BLIMP_stat_target```
 * Constraints and properties of the MPC can be modified in ```src/blimp_nmpc_wrapper_node/nodes/formation_config.py``` including the number of airships in the formation. Whenever this is changed, a recompilation is necessary.
 * ```submodules/AirCap/packages/3rdparty/airship_simulation/deflate_blimp.sh``` can be used to alter the airship buoyancy and rigidity
 * the script in ```submodules/AirCap/scripts/simulation/``` can be used to run various other experiments

### Real Blimp
To run the controller on a real blimp, the vehicle needs an OpenPilot Revolution Flight controller
and a companion computer with NVIDIA GPU connected via USB. We suggest a Jetson TX1 or TX2 with integrated wifi.
ROS should be run in a multi master setup using fkie_multimaster, with an instance of the MPC running on each vehicle.
