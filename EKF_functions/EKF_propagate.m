function [mxkm, Pxxkm] = EKF_propagate(dt, mxkm1, Pxxkm1, Pww, accel, gyro)

   g = [0; 0; 9.81]; % Gravity in NED frame
    
    % Extract state
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
    qxkm = qxkm1 + 0.5 * OM * qxkm1 * dt; % additive quat is a bad assumption
    qxkm = qxkm / norm(qxkm);  % Normalize quaternion

    accBkm = accBkm1;
    gyroBkm = gyroBkm1;
    
    mxkm = [pxkm; vxkm; qxkm; accBkm; gyroBkm];
   
    % -------- Jacobian F --------
    F = eye(16);
    
    % Position wrt velocity
    F(1:3, 4:6) = eye(3) * dt;
    
    % Velocity wrt quaternion
    F(4:6, 8:10) = -R * skew(a_corr) * dt;
    
    % Velocity wrt accel bias
    F(4:6, 11:13) = -R * dt; 
    
    % Quaternion wrt quaternion
    F(7:10, 7:10) = eye(4) + 0.5 * quatOmega(w_corr) * dt;
    
    % Quaternion wrt gyro bias
    Gq = quatVectorProductMatrix(qxkm1); 
    F(7:10, 14:16) = -0.5 * Gq * dt;

    % Covariance propagation
    Pxxkm = F * Pxxkm1 * F' + Pww;

end