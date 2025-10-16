# GNSS - IMU fusion

Extended Kalman filter (EKF) and Unscented Kalman filter (UKF) for loosely coupled GNSS - IMU fusion.

Next up: note on impact of process noise + comparison

Note: Processing the dataset all at once is slow. For practical problems, build using ROS. Example project: [GNSS-IMU-ErrorStateEKF](https://github.com/ram-bhaskara/GNSS-IMU-ErrorStateEKF) 

## EKF Results

Odometry from EKF (blue) vs truth (red)

![trajectory](/figures/EKF_trajectory.png)

EKF position estimates

![States](/figures/EKF_pos_states.png) 

## UKF Results

Odometry from UKF (blue) vs truth (red)

![trajectory](/figures/UKF_trajectory.png)

UKF position estimates

![States](/figures/UKF_pos_states.png) 

Dataset: [Zurich Urban Micro Aerial Vehicle Dataset](https://rpg.ifi.uzh.ch/zurichmavdataset.html)
