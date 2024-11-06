# 2D Pose Estimation for a Mobile Robot
Using an Extended Kalman Filter to reconstruct the 2D motion of a mobile robot

This was a course assignment for a State Estimation course at the Institute of Aerospace Studies at the University of Toronto. A robot drove around a lab environment with 17 landmarks equipped with reflective markers to generate a ground truth dataset through a motion capture system. An Extended Kalman Filter (EKF) was derived using the motion and observation models as well as histograms of the error in the sensors to recursively reconstruct the state of the robot at each timestep. The data used to reconstruct the states was the rotational and translational speed odometer measurements and laser rangefinder measurements with a 240° FOV of distance to the landmarks. 

The dataset used is property of UTIAS and is kept private, but in this repo I have attached my report for the EKF derivation and analysis of the results, an animation of the reconstructed motion, as well as all related MATLAB code.

![animation_gif](https://github.com/user-attachments/assets/1c86e763-7035-45d1-8d57-e686bcd04e6f)
