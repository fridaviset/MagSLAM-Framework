function [quat]=init_orient(u)

%Assume the sensor is stationary for a number of timesteps
%corresponding to the size of the inital window
initial_window=20;
quat=my_init_orientation(u,initial_window);

%Correct all initial orientation estimates with a factor that makes the 
%trajectory flat along the ground, and initially aligned with the x-axis.
%This rotation is found by tuning, and encoded in the variable eta_diff.
eta_diff=[-0.2; 0; 0.7];
quat=exp_q_L(eta_diff,quat);
quat=[quat(2:4);quat(1)];
end