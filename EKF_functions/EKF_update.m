function [mxkp, Pxxkp] = EKF_update(mxkm, Pxxkm, lla, v_meas, ref_lla, eph, epv, s_var, wgs84)
    
    ecef_ = lla2ecef([lla(1), lla(2), lla(3)], 'WGS84'); 
    [xN, yE, zD] = ecef2enu(ecef_(1), ecef_(2),ecef_(3), ...
                    ref_lla(1),ref_lla(2),ref_lla(3), ...
                    wgs84); 
    
    zk = [xN, yE, zD, v_meas]'; 
    mzkm = [mxkm(1:3); mxkm(4:6)]; 
    
    Hx = zeros(6, 16);
    Hx(1:3, 1:3) = eye(3);  
    Hx(4:6, 4:6) = eye(3);  
    
    Pvv = diag([eph^2, eph^2, epv^2, s_var^2, s_var^2, s_var^2]);
    
    Hv = 1; 
    
    Pxzkm = Pxxkm*Hx';
    Pzzkm = Hx*Pxxkm*Hx' + Hv*Pvv*Hv';
    Kk = Pxzkm/Pzzkm;
    mxkp = mxkm + Kk*(zk - mzkm);
    I = eye(size(Pxxkm));
    Pxxkp = (I - Kk * Hx) * Pxxkm * (I - Kk * Hx)' + Kk * Pvv * Kk';  % Joseph form
end