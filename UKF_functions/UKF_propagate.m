function [X_sigma_prop, mxkm, Pxxkm] = UKF_propagate(dt, mxkm1, Pxxkm1, Pww, accel, gyro, UKF_params)
    alpha = UKF_params(1);
    beta  = UKF_params(2);
    gamma_val = UKF_params(3);  
    lambda = UKF_params(4);
    n = numel(mxkm1);

    % sigma points
    X_sigma = X_sigma_points(mxkm1, Pxxkm1, gamma_val);
    L = size(X_sigma,2);
    
    % weights
    wm = [lambda/(n+lambda), repmat(0.5/(n+lambda),1,2*n)];
    wc = wm;
    wc(1) = wc(1) + (1 - alpha^2 + beta);
    
    % propagate each sigma through process model
    X_sigma_prop = zeros(n, L);
    for jj = 1:L
        X_sigma_prop(:,jj) = IMU_kinematics(dt, X_sigma(:,jj), accel, gyro);
    end
    
    % compute mean & cov
    mxkm = X_sigma_prop * wm';
    X1 = X_sigma_prop - mxkm(:, ones(1,L));
    Pxxkm = X1 * diag(wc) * X1' + Pww;
    
    Pxxkm = 0.5*(Pxxkm + Pxxkm');
end