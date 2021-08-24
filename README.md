# Integrated Model predictive control and Moving Horizon estimation for Stabilization at Goal location

In this repository, I demonstrate how to “Integrate Model predictive control and Moving Horizon estimation for Mobile robots” to achieve the two control objectives: Stabilization at Goal location.
My codes are further advanced implementation of the codes already explained in the workshop of M.W. Mehrez [1] because of following two reasons 1. I consider that both the estimates of current states and the applied control inputs are noisy which were not considered together in the examples of the workshop 2. I integrate together the Estimation and the Control problem which were addressed as two separate problems in the workshop. In fact it was also recommended in the workshop that “Combining MHE with MPC. A very good Exercise!” as the next step. I am grateful to M.W. Mehrez’s workshop which helped me to get started in the domain of solving Nonlinear Model predictive control using CASADI. 

# PROBLEM FORMULATION for Stabilization at Goal Location:
We have the below scenario for our first control objective i.e. Stabilization at the Goal location defined as x_s. The current position of the robot is not known and is measured by a range-bearing sensor which measures the ground truth states of the robot but has an unknown noise associated with it (assumed Gaussian here). The best estimate of the current position is determined by the MHE estimator which is then compared with the desired goal location for the MPC controller to determine the control trajectory for the entire pre-defined prediction horizon length. The first control action of the trajectory is applied to the low level controller of the robot which applies the higher level control action given to it by MPC but with an unknown noise associated with it (assumed Gaussian here). Therefore, the control input which ultimately goes to the robot is unknown, the effect of which can only be quantified by measuring states of the robot by range-bearing sensor which itself has noise associated with it. The goal is to make the robot reach the goal location under these Sensor and Actuator noise. The basic flowchart is given below for our problem. The important aspects of each of the Key blocks are explained and line no where they are implemented in the code MPC_MHE_Robot_PS_mul_shooting_vc1_stabilisation.m is also highlighted in the attached paper Paper_Integrated_MPC_MHE.pdf. 

![image](https://user-images.githubusercontent.com/83720464/130610808-939cd6e3-b6a6-4d00-8c72-c633b00f6665.png)
![image](https://user-images.githubusercontent.com/83720464/130616963-77e74a46-0d26-4143-8e4d-9d07b46abb64.png)
https://user-images.githubusercontent.com/83720464/130612103-a9ec6382-9521-4856-a3d3-279c908beff8.mp4

# REFERENCES
1.	M.W. Mehrez Optimization based solutions for control and State estimation in Dynamical systems (Implementation to Mobile Robots): A Workshop



