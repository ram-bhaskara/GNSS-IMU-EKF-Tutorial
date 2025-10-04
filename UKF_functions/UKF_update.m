function [mxkp, Pxxkp] = UKF_update(X_sigma_prop, mxkm, Pxxkm, zk, Pvv, UKF_params)
    alpha = UKF_params(1);
    beta  = UKF_params(2);
    lambda = UKF_params(4);
    
    n = numel(mxkm);
    L = size(X_sigma_prop,2);
    
    wm = [lambda/(n+lambda), repmat(0.5/(n+lambda),1,2*n)];
    wc = wm;
    wc(1) = wc(1) + (1 - alpha^2 + beta);

    % measurement sigma: here measurement is identity on first 6 states 
    Y_sigma = Y_sigma_points(X_sigma_prop);   % 6 x L
    mzkm = Y_sigma * wm';
    Y_diff = Y_sigma - mzkm(:,ones(1,L));
    Pzzkm = Y_diff * diag(wc) * Y_diff' + Pvv;
    x_diff = X_sigma_prop - mxkm(:,ones(1,L));
    Pxzkm = x_diff * diag(wc) * Y_diff';
    Kk = Pxzkm * pinv(Pzzkm);
    mxkp = mxkm + Kk * (zk - mzkm);
    Pxxkp = Pxxkm - Pxzkm * Kk' - Kk * Pxzkm' + Kk * Pzzkm * Kk';
    Pxxkp = 0.5*(Pxxkp + Pxxkp');
end