function mxkm = IMU_kinematics(dt, mxkm1, accel, gyro)
    
    mxkm1 = reshape(mxkm1, [numel(mxkm1), 1]); % state vector 
    g = [0; 0; -9.81]; % Gravity in NED frame

    pxkm1     = mxkm1(1:3); 
    vxkm1     = mxkm1(4:6); 
    qxkm1     = mxkm1(7:10);  % quaternion (scalar-first)
    accBkm1   = mxkm1(11:13); 
    gyroBkm1  = mxkm1(14:16); 
    
    % Sensor correction
    a_corr = accel - accBkm1;
    w_corr = gyro  - gyroBkm1;
    
    % Rotation matrix (body to world)
    R = quat2rotm(qxkm1');  % transpose needed for MATLAB

    % Propagate state
    pxkm = pxkm1 + vxkm1 * dt; 
    vxkm = vxkm1 + (R * a_corr + g) * dt;

    OM = quatOmega(w_corr); 
    qxkm = qxkm1 + 0.5 * OM * qxkm1 * dt;
    qxkm = qxkm / norm(qxkm);  % Normalize quaternion

    accBkm = accBkm1;
    gyroBkm = gyroBkm1;
    
    mxkm = [pxkm; vxkm; qxkm; accBkm; gyroBkm];
    
end