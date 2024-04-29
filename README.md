# MATLAB_ICP
Various ICP algorithm implementation with MATLAB

1. Click /point_cloud/PointCloudwithNoise.mat to load point cloud data
2. Set Configure

```
lambda = 0.0; % lambda for LM optimization. If 0, Gauss-Newton
max_iter = 50;
trans_epsilon = 0.001; % meter
rot_epsilone = 0.1; % deg
max_distance = 0.2; % Not used now

ndt_grid_size = 0.2; % Not used now

debug = true; % If true, show ICP progress in visualiation

% Target, Moving ground truth transformation
gt_trans_x = 1.0;
gt_trans_y = 1.5;
gt_rot_yaw = 25 * pi/180; 
```

