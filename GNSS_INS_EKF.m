clc; 
% clear; 

dataPath = '.\'; 
addpath(dataPath); 
%% DATA READ
[accelData, gyroData, GPSData] = dataRead(); 
%% Init

% state: [r, v, q, b_a, b_g]
% initial position from GPS fix

% wgs84 = wgs84Ellipsoid('meter'); 
ref_lla = GPSData(1,2:4);
% ecef_r0 = lla2ecef(GPSData(1,2:4));
% [x0,y0,z0] = ecef2ned(ecef2(1), ecef2(2), ecef2(3), GPSData(1,2), GPSData(1,3), GPSData(1,4), wgs84); 
% r0 = [x0, y0, z0];
[posx, posy, posz] = gps2cart(ref_lla(1), ref_lla(2), ref_lla(3)); 
r0 = [posx, posy, posz]; 
% r0 = ecef_r0; 
v0_ned = GPSData(1,5:7);
% [v0x, v0y, v0z] = ned2ecef(v0_ned(1), v0_ned(2),v0_ned(3), ref_lla(1), ref_lla(2), ref_lla(3), 'wgs84');
% v0 = [v0x, v0y, v0z]; 
v0 = v0_ned;

q0 = [1 0 0 0]; 
accelBias0 = [0 0 0]; 
gyroBias0 = [0 0 0];


mx0 = [r0 v0 q0 accelBias0 gyroBias0]';
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

%% Propagation step - 
% Asynchronous accel and gyro data
imuTimes = unique([accelData(:,1)]);
gpsTimes = GPSData(:,1);

allTimes = sort(unique([imuTimes; gpsTimes]));

% EKF

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
    
    [mxkm, Pxxkm] = EKF_propagate(dt, mxkm1, Pxxkm1, Pww, acc_k, gyro_k); 

    gpsIdx = find(gpsTimes == tk, 1);
    if ~isempty(gpsIdx)
        lla = GPSData(gpsIdx, 2:4);             % [lat, lon, alt]
        v_meas = GPSData(gpsIdx, 5:7);         % [vn; ve; vd]
        eph  = GPSData(gpsIdx, 8);              % horizontal DOP
        epv  = GPSData(gpsIdx, 9);              % vertical DOP
        s_var = GPSData(gpsIdx, 10);            % velocity variance

        [mxkp, Pxxkp] = EKF_update(mxkm, Pxxkm, lla, v_meas, ref_lla, eph, epv, s_var);
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
figure; 
x_est = mxstore(1,:); 
y_est = mxstore(2,:); 
z_est = mxstore(3,:); 

vx_est = mxstore(4,:); 
vy_est = mxstore(5,:); 
vz_est = mxstore(6,:); 

txstore = (txstore - txstore(1))*1e-6; 
subplot(3,1,1)
plot(txstore, x_est); ylabel('x')
% hold on; 
% plot(txstore, x_est+sxstore(1,:), 'r', txstore, -sxstore(1,:), 'r'); hold off
subplot(3,1,2)
plot(txstore, y_est); ylabel('y')
% hold on; 
% plot(txstore, sxstore(2,:), 'r', txstore, -sxstore(2,:), 'r'); hold off
subplot(3,1,3)
plot(txstore, z_est); ylabel('z')
% hold on; 
% plot(txstore, sxstore(3,:), 'r', txstore, -sxstore(3,:), 'r'); hold off

figure; 
subplot(3,1,1)
plot(x_gps); ylabel('x')
subplot(3,1,2)
plot(y_gps); ylabel('y')
subplot(3,1,3)
plot(z_gps); ylabel('z')
title("Onboard GPS")

%
figure; 
% plot3(x_est, y_est, z_est, 'r'); 
% hold on; 
% plot3(x_gt, y_gt, z_gt, 'b'); 
% hold off;
%% Functions

function OM = quatOmega(omega)
    wx = omega(1); wy = omega(2); wz = omega(3);
    OM = [  0   -wx -wy -wz;
           wx    0   wz -wy;
           wy  -wz   0   wx;
           wz   wy -wx   0];
end

function Gq = G(q)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    Gq = [ -q1, -q2, -q3;   % Row 1
            q0, -q3,  q2;   % Row 2
            q3,  q0, -q1;   % Row 3
           -q2,  q1,  q0];  % Row 4
end

function S = skew(v)
    S = [   0   -v(3)  v(2);
          v(3)   0   -v(1);
         -v(2)  v(1)   0];
end

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
    qxkm = qxkm1 + 0.5 * OM * qxkm1 * dt;
    qxkm = qxkm / norm(qxkm);  % Normalize quaternion

    accBkm = accBkm1;
    gyroBkm = gyroBkm1;
    
    mxkm = [pxkm; vxkm; qxkm; accBkm; gyroBkm];
   
    % -------- Jacobian F --------
    F = eye(16);
    
    % Position wrt velocity
    F(1:3, 4:6) = eye(3) * dt;
    
    % Velocity wrt quaternion (∂v/∂q ≈ -R * skew(a_corr) * dt)
    F(4:6, 8:10) = -R * skew(a_corr) * dt;
    
    % Velocity wrt accel bias
    F(4:6, 11:13) = -R * dt; 
    
    % Quaternion wrt quaternion
    F(7:10, 7:10) = eye(4) + 0.5 * quatOmega(w_corr) * dt;
    
    % Quaternion wrt gyro bias
    Gq = G(qxkm1); 
    F(7:10, 14:16) = -0.5 * Gq * dt;

    % Covariance propagation
    Pxxkm = F * Pxxkm1 * F' + Pww;

end

function [mxkp, Pxxkp] = EKF_update(mxkm, Pxxkm, lla, v_meas, ref_lla, eph, epv, s_var)
    
ecef = lla2ecef(lla);         % MATLAB function
wgs84 = wgs84Ellipsoid('meter');
% [px, py, pz] = ecef2ned(ecef(1), ecef(2), ecef(3), ...
%                 ref_lla(1), ref_lla(2), ref_lla(3), wgs84); 

% p_meas  = [px; py; pz];

[posx, posy, posz] = gps2cart(lla(1), lla(2), lla(3)); 
ecef = [posx, posy, posz]; 
% [v0x, v0y, v0z] = ned2ecef(v_meas(1), v_meas(2),v_meas(3), ref_lla(1), ref_lla(2), ref_lla(3), wgs84);
% v_ecef = [v0x, v0y, v0z];
% zk = [p_meas; v_meas]; % in ned frame
zk = [ecef, v_meas]'; 
mzkm = [mxkm(1:3); mxkm(4:6)]; 

Hx = zeros(6, 16);
Hx(1:3, 1:3) = eye(3);  % 
Hx(4:6, 4:6) = eye(3);  % 

Pvv = diag([eph^2, eph^2, epv^2, s_var^2, s_var^2, s_var^2]);

Hv = 1; 

    
    Pxzkm = Pxxkm*Hx';
    Pzzkm = Hx*Pxxkm*Hx' + Hv*Pvv*Hv';
    Kk = Pxzkm/Pzzkm;
    mxkp = mxkm + Kk*(zk - mzkm);
    I = eye(size(Pxxkm));
    Pxxkp = (I - Kk * Hx) * Pxxkm * (I - Kk * Hx)' + Kk * Pvv * Kk';  % Joseph form
end

function [x,y,z] = gps2cart(lat,lon,alt)

    % WGS84 ellipsoid constants:
    a = 6378137; % earth radius (semi-major axis of the ellipse)
    e = 8.1819190842622e-2; % first eccentricity
    
    % intermediate calculation
    % N is vertical radius of curvature
    
    lat = lat / 180 * pi;
    lon = lon / 180 * pi;
    N = a ./ sqrt(1 - e^2 .* sin(lat).^2);
    
    % results:
%     x = (N + alt) .* cos(lat) .* cos(lon);
    x = (N + alt) .* cos(lat) .* lon;
%     y = (N + alt) .* cos(lat) .* sin(lon);
    y = (N + alt) .* lat;
%     z = ((1 - e^2) .* N + alt) .* sin(lat);
    z = alt;
end
