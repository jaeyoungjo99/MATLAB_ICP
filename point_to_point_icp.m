%%  File Name: point_to_point_icp.m %%%%%%%%%%%%%%%%%%%%
%
%  Description: Point to Point ICP code level implementation
%
%  By Jaeyoung Jo, AI LAB (wodud3743@gmail.com) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Loading Point Cloud & Configure

% Configure
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

% Loading Point Cloud
data = load('PointCloudwithNoise.mat');  % 'example.mat' 파일에서 데이터 불러오기
pt_cloud_origin = data.inputCloud;
pt_cloud_origin = pcdownsample(pt_cloud_origin,'gridAverage',0.01);

method = "point2plane"; % svd, point2point, point2plane, ndt, gicp

% affine3d input matrix is transpose of general R matrix
gt_transform_matrix = [cos(gt_rot_yaw) sin(gt_rot_yaw) 0 0; ...
                 -sin(gt_rot_yaw) cos(gt_rot_yaw) 0 0; ...
                 0 0 1 0; ...
                 gt_trans_x gt_trans_y 0 1];

gt_transform_affine = affine3d(gt_transform_matrix);

% Transform moving point cloud
pt_cloud_moving = pctransform(pt_cloud_origin,gt_transform_affine);

% Downsampling
pt_cloud_moving_ds = pcdownsample(pt_cloud_moving,'gridAverage',0.02);

%% Original Display

% Display the Original and Affine Transformed 3-D Point Clouds
figure1 = figure('WindowState','normal');
axes1 = axes('Parent',figure1);
pcshow(pt_cloud_origin,'Parent',axes1); 
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3-D Point Cloud','FontSize',14)

pcshowpair(pt_cloud_origin,pt_cloud_moving_ds,'VerticalAxis','Z','VerticalAxisDir','Up')


%% Run Registration
fprintf('Target point Num %d \n',pt_cloud_origin.Count);
fprintf('Moving point Num %d \n',pt_cloud_moving_ds.Count);

pt_cloud_transformed = pt_cloud_moving_ds;

% SVD
if method == "svd"
    % 1. SVD Based
    for iter = 1:max_iter
        fprintf("Iter %d \n",iter);
        
        % Finding corresponding
        numPoints = pt_cloud_transformed.Count;
        nearestIndices = zeros(numPoints, 1);
        
        for i = 1:numPoints
            % pc1의 i번째 점에 대해 pc2에서 가장 가까운 한 개의 점 찾기
            [indices, distances] = findNearestNeighbors(pt_cloud_origin, pt_cloud_transformed.Location(i,:), 1);
            
            % 가장 가까운 점의 인덱스와 거리 저장
            nearestIndices(i) = indices;
        end
    
        pt_cloud_matched = select(pt_cloud_origin, nearestIndices);
        fprintf('Target point matched Num %d \n',pt_cloud_matched.Count);
        
        H = zeros(3,3);
        
        target_mean = [mean(pt_cloud_matched.Location(:,1));
                        mean(pt_cloud_matched.Location(:,2));
                        mean(pt_cloud_matched.Location(:,3));];
        
        moving_mean = [mean(pt_cloud_transformed.Location(:,1));
                        mean(pt_cloud_transformed.Location(:,2));
                        mean(pt_cloud_transformed.Location(:,3));];
        
        
        for i = 1:numPoints
            target_displaced = pt_cloud_matched.Location(i,:)' - target_mean;
            moving_displaced = pt_cloud_transformed.Location(i,:)' - moving_mean;
        
            H = H + moving_displaced*target_displaced';
        end
        
        [U, S, V] = svd(H);
        
        R = V'*U;
        
        pitch = asin(-R(3,1));
        if cos(pitch) ~= 0
            roll = atan2(R(3,2)/cos(pitch), R(3,3)/cos(pitch));
            yaw = atan2(R(2,1)/cos(pitch), R(1,1)/cos(pitch));
        else
            roll = atan2(-R(1,2), R(1,3));
            yaw = 0;  % 자유도 상실로 정확한 값 계산 불가
        end
        
        % 각도를 도(degree)로 변환 (선택적)
        roll = rad2deg(roll);
        pitch = rad2deg(pitch);
        yaw = rad2deg(yaw);
        
        T = target_mean - R*moving_mean;
        
        fprintf('Roll %f, Pitch %f, Yaw %f \n',roll, pitch, yaw);
        fprintf('X %f, Y %f, Z %f \n',T(1), T(2), T(3));
    
        pt_transform = [R(1,1) R(1,2) R(1,3) 0; ...
                         R(2,1) R(2,2) R(2,3) 0; ...
                         R(3,1) R(3,2) R(3,3) 0; ...
                         T(1) T(2) T(3) 1];
    
        tform_iter = affine3d(pt_transform);
        
        % Transform moving point cloud
        pt_cloud_transformed = pctransform(pt_cloud_transformed,tform_iter);
    
    end

% Point to Point
elseif method == "point2point"
    total_transform = zeros(3,1);
    for iter = 1:max_iter
        fprintf("P2Point Iter %d \n",iter);
        
        % Finding corresponding
        numPoints = pt_cloud_transformed.Count;
        coeffs = zeros(numPoints,3);
        coeff_distances = zeros(numPoints,1);

        
        
        for i = 1:numPoints
            % pc1의 i번째 점에 대해 pc2에서 가장 가까운 한 개의 점 찾기
            [indices, distances] = findNearestNeighbors(pt_cloud_origin, pt_cloud_transformed.Location(i,:), 1);
            
            % 가장 가까운 점의 인덱스와 거리 저장
            coeff_distances(i) = distances;

            coeffs(i,:) = pt_cloud_transformed.Location(i,:)' - pt_cloud_origin.Location(indices,:)';
            coeff_norm = sqrt(coeffs(i,1)*coeffs(i,1) + coeffs(i,2)*coeffs(i,2));

            coeffs(i,1) = coeffs(i,1) / coeff_norm;
            coeffs(i,2) = coeffs(i,2) / coeff_norm;
        end


        matA = zeros(numPoints,3);
        matAt = zeros(3,numPoints);
        matAtA = zeros(3,3);
        matB = zeros(numPoints,1);
        matAtB = zeros(3,1);
        matAtAdiag = zeros(3,3);
        matX = zeros(3,1); % yaw x y 
      
        % Calculate residual vector
        for i = 1:numPoints
            rot_yaw = pt_cloud_transformed.Location(i,1) * (coeffs(i,2)) ...
                    - pt_cloud_transformed.Location(i,2) * (coeffs(i,1));

            matA(i,1) = rot_yaw;
            matA(i,2) = coeffs(i,1);
            matA(i,3) = coeffs(i,2);
            matB(i,1) = -coeff_distances(i);
        end

        matAt = matA';
        
        diagonal = diag(matAt * matA);
        matAtAdiag = diag(diagonal);

        matAtA = matAt * matA + lambda * matAtAdiag;
        matAtB = matAt * matB; 
        
        matX = matAtA \ matAtB;

        fprintf("Dyaw: %f, dx %f, dy %f \n",matX(1)*180/pi, matX(2), matX(3));

        pt_transform = [cos(matX(1)) sin(matX(1)) 0 0; ...
                       -sin(matX(1)) cos(matX(1)) 0 0; ...
                         0 0 1 0; ...
                         matX(2) matX(3) 0 1];
    
        tform_iter = affine3d(pt_transform);
        
        % Transform moving point cloud
        pt_cloud_transformed = pctransform(pt_cloud_transformed,tform_iter);

        total_transform = total_transform + matX;

        if abs(matX(2)) < trans_epsilon && abs(matX(3)) < trans_epsilon && abs(matX(1))*180/pi < rot_epsilone
            break
        end

        if debug == true
            pcshowpair(pt_cloud_origin,pt_cloud_transformed,'VerticalAxis','Z','VerticalAxisDir','Up');
            pause(0.1);
        end

    end
    fprintf("Total Dyaw: %f, dx %f, dy %f \n",total_transform(1)*180/pi, total_transform(2), total_transform(3));
    fprintf("Error Dyaw: %f, dx %f, dy %f \n",(gt_rot_yaw + total_transform(1))*180/pi, gt_trans_x + total_transform(2), gt_trans_y + total_transform(3));
  
% Point to Plane
elseif method == "point2plane"
    total_transform = zeros(3,1);
    affine_matrix = [1 0 0 0; ...
                    0 1 0 0; ...
                    0 0 1 0; ...
                    0 0 0 1];
    for iter = 1:max_iter
        fprintf("P2Plane Iter %d \n",iter);
        
        % Finding corresponding
        numPoints = pt_cloud_transformed.Count;
        coeffs = zeros(numPoints,3);
        coeff_distances = zeros(numPoints,1);
        
        for i = 1:numPoints
            % pc1의 i번째 점에 대해 pc2에서 가장 가까운 한 개의 점 찾기
            [indices, distances] = findNearestNeighbors(pt_cloud_origin, pt_cloud_transformed.Location(i,:), 5);
            % 가장 가까운 점의 인덱스와 거리 저장

            matA0 = zeros(5, 3);
            % matB0 = zeros(5, 1)
            matB0 = [-1; -1; -1; -1; -1];
            for j = 1:5
                matA0(j, 1) = pt_cloud_origin.Location(indices(j),1);
                matA0(j, 2) = pt_cloud_origin.Location(indices(j),2);
                matA0(j, 3) = pt_cloud_origin.Location(indices(j),3);
            end
            matX0 = matA0 \ matB0;
            pa = matX0(1);  % x계수
            pb = matX0(2);  % y계수
            pc = matX0(3);  % z계수
            pd = 1;
            ps = sqrt(pa * pa + pb * pb + pc * pc);
            pa = pa/ps; pb = pb/ps; pc = pc/ps; pd = pd/ps;

            pd2 = pa * pt_cloud_transformed.Location(i,1) ...
                + pb * pt_cloud_transformed.Location(i,2) ...
                + pc * pt_cloud_transformed.Location(i,3) + pd;
            weight = 20;
            coeffs(i,:) = [pa pb pc] * weight;
            coeff_distances(i) = pd2 * weight;
        end

        matA = zeros(numPoints,3);
        matAt = zeros(3,numPoints);
        matAtA = zeros(3,3);
        matB = zeros(numPoints,1);
        matAtB = zeros(3,1);
        matAtAdiag = zeros(3,3);
        matX = zeros(3,1); % yaw x y 
      
        % Calculate residual vector
        for i = 1:numPoints

            rot_yaw = pt_cloud_transformed.Location(i,1) * (coeffs(i,2)) ...
                    - pt_cloud_transformed.Location(i,2) * (coeffs(i,1));

            matA(i,1) = rot_yaw;
            matA(i,2) = coeffs(i,1);
            matA(i,3) = coeffs(i,2);
            matB(i,1) = -coeff_distances(i);
        end

        % LM Optimization
        matAt = matA';
        
        diagonal = diag(matAt * matA);
        matAtAdiag = diag(diagonal);

        matAtA = matAt * matA + lambda * matAtAdiag;
        matAtB = matAt * matB; 
        
        matX = matAtA \ matAtB;

        fprintf("Dyaw: %f, dx %f, dy %f \n",matX(1)*180/pi, matX(2), matX(3));

        pt_transform = [cos(matX(1)) sin(matX(1)) 0 0; ...
                       -sin(matX(1)) cos(matX(1)) 0 0; ...
                         0 0 1 0; ...
                         matX(2) matX(3) 0 1];
    
        tform_iter = affine3d(pt_transform);

        affine_matrix = affine_matrix * pt_transform;
        
        % Transform moving point cloud
        pt_cloud_transformed = pctransform(pt_cloud_transformed,tform_iter);

        total_transform = total_transform + matX;

        if abs(matX(2)) < trans_epsilon && abs(matX(3)) < trans_epsilon && abs(matX(1))*180/pi < rot_epsilone
            break
        end

        if debug == true
            pcshowpair(pt_cloud_origin,pt_cloud_transformed,'VerticalAxis','Z','VerticalAxisDir','Up');
            pause(0.1);
        end
    end
    
    translation = affine_matrix(4, 1:3);
    rotationMatrix = affine_matrix(1:3, 1:3);
    end_pitch = asin(-rotationMatrix(3,1));
    if cos(end_pitch) ~= 0
        end_roll = atan2(rotationMatrix(3,2)/cos(end_pitch), rotationMatrix(3,3)/cos(end_pitch));
        end_yaw = atan2(rotationMatrix(2,1)/cos(end_pitch), rotationMatrix(1,1)/cos(end_pitch));
    else
        end_roll = atan2(-rotationMatrix(1,2), rotationMatrix(1,3));
        end_yaw = 0;  % 자유도 상실
    end

    fprintf("Total Dyaw: %f, dx %f, dy %f \n",end_yaw*180/pi, translation(1), translation(2));
    fprintf("Error Dyaw: %f, dx %f, dy %f \n",(gt_rot_yaw + end_yaw)*180/pi, gt_trans_x + translation(1), gt_trans_y + translation(2));

% NDT (TODO)
elseif method == "ndt"
    total_transform = zeros(3,1);
    
    % Map Generation

    % Iteration start
    for iter = 1:max_iter
    
    end

    fprintf("Total Dyaw: %f, dx %f, dy %f \n",total_transform(1)*180/pi, total_transform(2), total_transform(3));
    fprintf("Error Dyaw: %f, dx %f, dy %f \n",(gt_rot_yaw + total_transform(1))*180/pi, gt_trans_x + total_transform(2), gt_trans_y + total_transform(3));

% GICP
elseif method == "gicp"
    total_transform = zeros(3,1);
    affine_matrix = [1 0 0 0; ...
                    0 1 0 0; ...
                    0 0 1 0; ...
                    0 0 0 1];
    sourceNumPoints = pt_cloud_transformed.Count;
    tartgetNumPoints = pt_cloud_origin.Count;

    source_covariances = zeros(3,3,sourceNumPoints);
    target_covariances = zeros(3,3,tartgetNumPoints);

    source_cov_rotations = zeros(3,3,sourceNumPoints);
    target_cov_rotations = zeros(3,3,tartgetNumPoints);

    weights = zeros(sourceNumPoints,1);
    
    e = 0.3;
    C = [1 0 0;
         0 1 0;
         0 0 e];

    % Target point covariance in advance
    for i = 1:tartgetNumPoints
        [indices, distances] = findNearestNeighbors(pt_cloud_origin, pt_cloud_origin.Location(i,:), 5);
        neighbors = pt_cloud_origin.Location(indices, :);
        meanPoint = mean(neighbors, 1);
        centeredPoints = neighbors - meanPoint;
        covarianceMatrix = (centeredPoints' * centeredPoints) / (size(neighbors, 1) - 1);
        [V, D] = eig(covarianceMatrix);
        target_cov_rotations(:,:,i) = V;
    end

    for iter = 1:max_iter
        fprintf("GICP Iter %d \n",iter);
        
        % Finding corresponding
        
        coeffs = zeros(sourceNumPoints,3);
        coeff_distances = zeros(sourceNumPoints,1);
        
        for i = 1:sourceNumPoints
            % source point covariance
            [source_indices, source_distances] = findNearestNeighbors(pt_cloud_transformed, pt_cloud_transformed.Location(i,:), 5);
            neighbors = pt_cloud_transformed.Location(source_indices, :);
            meanPoint = mean(neighbors, 1);
            centeredPoints = neighbors - meanPoint;
            covarianceMatrix = (centeredPoints' * centeredPoints) / (size(neighbors, 1) - 1);
            [V, D] = eig(covarianceMatrix);
            source_cov_rotations(:,:,i) = V;
            
            % find nearest target covariance
            [target_indice, target_distance] = findNearestNeighbors(pt_cloud_origin, pt_cloud_transformed.Location(i,:), 1);
            
            source_rotation = source_cov_rotations(:,:,i);
            target_rotation = target_cov_rotations(:,:,target_indice);

            
            
          
            dist_vec = (pt_cloud_transformed.Location(i,:)' - pt_cloud_origin.Location(target_indice,:)');
                        coeff_norm = sqrt(dist_vec(1)*dist_vec(1) + dist_vec(2)*dist_vec(2));

            Ca = source_rotation*C*source_rotation';
            Cb = target_rotation*C*target_rotation';
            
            weight = dist_vec'*inv(Ca + Cb) * dist_vec;
            weights(i) = weight;
            
            coeffs(i,:) = dist_vec / coeff_norm * weight;
            coeff_distances(i) = target_distance * weight;
        end

        matA = zeros(sourceNumPoints,3);
        matAt = zeros(3,sourceNumPoints);
        matAtA = zeros(3,3);
        matB = zeros(sourceNumPoints,1);
        matAtB = zeros(3,1);
        matAtAdiag = zeros(3,3);
        matX = zeros(3,1); % yaw x y 
      
        % Calculate residual vector
        for i = 1:sourceNumPoints

            rot_yaw = pt_cloud_transformed.Location(i,1) * (coeffs(i,2)) ...
                    - pt_cloud_transformed.Location(i,2) * (coeffs(i,1));

            matA(i,1) = rot_yaw;
            matA(i,2) = coeffs(i,1);
            matA(i,3) = coeffs(i,2);
            matB(i,1) = -coeff_distances(i);
        end

        % LM Optimization
        matAt = matA';
        
        diagonal = diag(matAt * matA);
        matAtAdiag = diag(diagonal);

        matAtA = matAt * matA + lambda * matAtAdiag;
        matAtB = matAt * matB; 
        
        matX = matAtA \ matAtB;

        fprintf("Dyaw: %f, dx %f, dy %f \n",matX(1)*180/pi, matX(2), matX(3));

        pt_transform = [cos(matX(1)) sin(matX(1)) 0 0; ...
                       -sin(matX(1)) cos(matX(1)) 0 0; ...
                         0 0 1 0; ...
                         matX(2) matX(3) 0 1];
    
        tform_iter = affine3d(pt_transform);

        affine_matrix = affine_matrix * pt_transform;
        
        % Transform moving point cloud
        pt_cloud_transformed = pctransform(pt_cloud_transformed,tform_iter);

        total_transform = total_transform + matX;

        if abs(matX(2)) < trans_epsilon && abs(matX(3)) < trans_epsilon && abs(matX(1))*180/pi < rot_epsilone
            break
        end

        if debug == true
            pcshowpair(pt_cloud_origin,pt_cloud_transformed,'VerticalAxis','Z','VerticalAxisDir','Up');
            pause(0.1);
        end
    end
    
    translation = affine_matrix(4, 1:3);
    rotationMatrix = affine_matrix(1:3, 1:3);
    end_pitch = asin(-rotationMatrix(3,1));
    if cos(end_pitch) ~= 0
        end_roll = atan2(rotationMatrix(3,2)/cos(end_pitch), rotationMatrix(3,3)/cos(end_pitch));
        end_yaw = atan2(rotationMatrix(2,1)/cos(end_pitch), rotationMatrix(1,1)/cos(end_pitch));
    else
        end_roll = atan2(-rotationMatrix(1,2), rotationMatrix(1,3));
        end_yaw = 0;  % 자유도 상실
    end

    fprintf("Total Dyaw: %f, dx %f, dy %f \n",end_yaw*180/pi, translation(1), translation(2));
    fprintf("Error Dyaw: %f, dx %f, dy %f \n",(gt_rot_yaw + end_yaw)*180/pi, gt_trans_x + translation(1), gt_trans_y + translation(2));
end

%% Display the Original and Affine Transformed 3-D Point Clouds

pcshowpair(pt_cloud_origin,pt_cloud_transformed,'VerticalAxis','Z','VerticalAxisDir','Up')


