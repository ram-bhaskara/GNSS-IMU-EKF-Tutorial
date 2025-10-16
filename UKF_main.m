clc; clear; 
%% DATA READ
addpath('./UKF_functions/','math_utils');
dataPath = './data';
[accelData, gyroData, GPSData] = dataRead(dataPath); 
run('loadGroundTruthAGL.m'); % ground truth - their definition
%%
ref_lla = GPSData(1,2:4);
pos_ecef = lla2ecef([ref_lla(1), ref_lla(2), ref_lla(3)], 'WGS84');
r0 = pos_ecef; 

% ecef2enu
wgs84 = wgs84Ellipsoid('meter');
[xE, yN, zU] = ecef2enu(r0(1),r0(2),r0(3),ref_lla(1),ref_lla(2),ref_lla(3),wgs84); 

v0_ned = GPSData(1,5:7);
v0 = v0_ned;

q0 = [1 0 0 0]; 
accelBias0 = [0 0 0]; 
gyroBias0 = [0 0 0];

g = [0; 0; 9.81];

mx0 = [zeros(1,3) v0 q0 accelBias0 gyroBias0]';
nStates = length(mx0); 

eph = GPSData(1,8); epv = GPSData(1,9); %hDop, vDop in meters
s_var = GPSData(1,10); 

Prr0 = diag([eph^2, eph^2, epv^2]); 
Pvv0 = s_var^2 * eye(3); 
Pqq0 = 	0.01^2 * eye(4); 
Pab0 = 1e-3 * eye(3); 
Pgb0 = 1e-4 * eye(3);

Pxx0 = blkdiag(Prr0, Pvv0, Pqq0, Pab0, Pgb0);
Pww = Pxx0; % changes based on accel readings

Q_pos = 1e-3 * eye(3);
Q_vel = 1e-3 * eye(3);
Q_q   = 1e-4 * eye(4);
Q_ab  = 1e-5 * eye(3);
Q_gb  = 1e-5 * eye(3);
Pww = blkdiag(Q_pos, Q_vel, Q_q, Q_ab, Q_gb);

%% UKF Parameters
% UKF parameters 
alpha = 0.5; 
k = 0; 
beta = 2;
n = numel(mx0); % numel(Hx) = n (State dimension)
gamma = alpha^2 * (n + k); 
lambda = gamma - n;
gamma = sqrt(gamma); 
UKF_params = [alpha, beta, gamma, lambda]; 
%% Propagation
% Asynchronous accel and gyro data
imuTimes = unique([accelData(:,1)]);
gpsTimes = GPSData(:,1);

allTimes = sort(unique([imuTimes; gpsTimes]));

mxkm1 = mx0; 
Pxxkm1 = Pxx0;

xcount = 1;
txstore(:,xcount) = allTimes(1);
mxstore(:,xcount) = mxkm1;
sxstore(:,xcount) = sqrt(diag(Pxxkm1));

for k = 2:length(allTimes)
    tk = allTimes(k); 
    tk1 = allTimes(k-1); 
    dt = (tk - tk1) * 1e-6; 
    
    % --- Get closest IMU data (or interpolate if you prefer)
    [~, idxA] = min(abs(accelData(:,1) - tk));
    [~, idxG] = min(abs(gyroData(:,1)  - tk));
    acc_k = accelData(idxA, 2:4)';  % [ax; ay; az]
    gyro_k = gyroData(idxG, 2:4)';  % [gx; gy; gz]
    
    [X_sigma_prop, mxkm, Pxxkm] = UKF_propagate(dt, mxkm1, Pxxkm1, Pww, acc_k, gyro_k, UKF_params); 
    Pxxkm = 0.5*(Pxxkm + Pxxkm');

    gpsIdx = find(gpsTimes == tk, 1);
    if ~isempty(gpsIdx)
        lla = GPSData(gpsIdx, 2:4);             % [lat, lon, alt]
        v_meas = GPSData(gpsIdx, 5:7);         % [vn; ve; vu]
        eph  = GPSData(gpsIdx, 8);              % horizontal DOP
        epv  = GPSData(gpsIdx, 9);              % vertical DOP
        s_var = GPSData(gpsIdx, 10);            % velocity variance

        Pvv = diag([eph^2, eph^2, epv^2, s_var^2, s_var^2, s_var^2]);
       
        ecef_ = lla2ecef([lla(1), lla(2), lla(3)], 'WGS84'); 
        [xN, yE, zD] = ecef2enu(ecef_(1), ecef_(2),ecef_(3), ...
                        ref_lla(1),ref_lla(2),ref_lla(3), ...
                        wgs84); 
        
        zk = [xN, yE, zD, v_meas]'; 
        [mxkp, Pxxkp] = UKF_update(X_sigma_prop, mxkm, Pxxkm, zk, Pvv, UKF_params);
    else
        mxkp = mxkm;
        Pxxkp = Pxxkm;
    end

    xcount = xcount + 1;
    txstore(xcount) = tk;
    mxstore(:,xcount) = mxkp;
    sxstore(:,xcount) = sqrt(diag(Pxxkp));
    
    % === Recycle for next iteration
    mxkm1 = mxkp;
    Pxxkm1 = Pxxkp;
end

%% Plots
x_est = mxstore(1,:); 
y_est = mxstore(2,:); 
z_est = mxstore(3,:); 
vx_est = mxstore(4,:); 
vy_est = mxstore(5,:); 
vz_est = mxstore(6,:); 

tx = (txstore - txstore(1))*1e-6; 

x_GT = x_gps - x_gps(1);
y_GT = y_gps - y_gps(1); 
z_GT = z_gps - z_gps(1); 

figure; 
subplot(3,1,1)
plot(tx, x_est,'lineWidth', 2.5); ylabel('x [m]');
hold on; 
plot(x_GT, '-.r','lineWidth',1.5); 
grid on; axis tight; legend('UKF_x','Truth')

subplot(3,1,2)
plot(tx, y_est,'lineWidth', 2.5); ylabel('y [m]')
hold on; 
plot(y_GT, '-.r','lineWidth',1.5);
grid on; axis tight; legend('UKF_y','Truth')

subplot(3,1,3)
plot(tx, z_est,'lineWidth', 2.5); ylabel('z [m]')
hold on; 
plot(z_GT, '-.r','lineWidth',1.5); 
grid on; axis tight; legend('UKF_z','Truth')

%
figure; 
plot3(x_est, y_est, z_est,'lineWidth',2.5); 
hold on; 
plot3(x_GT, y_GT, z_GT, '-.r','lineWidth',3); 
hold off;
title('UKF Trajectory: Estimate vs Ground Truth'); 
xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]'); 
grid on; 
legend('Estimated position', 'True trajectory'); 
view([0 90]);