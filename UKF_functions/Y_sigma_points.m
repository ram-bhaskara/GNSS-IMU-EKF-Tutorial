function Y = Y_sigma_points(X_sigma_prop)
    % measurement is first 6 state components: [x; y; z; vx; vy; vz] in ENU
    Y = X_sigma_prop(1:6, :);
end