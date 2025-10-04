function X = X_sigma_points(mxkm1, Pxxkm1, gamma_val)
    n = numel(mxkm1);
    S = chol(Pxxkm1, 'lower');   % lower triangular L s.t. L*L' = P
    Y = mxkm1(:, ones(1,n));
    X = [mxkm1, Y + gamma_val * S, Y - gamma_val * S];
end