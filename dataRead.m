function [accelData, gyroData, GPSData] = dataRead()

accelFile = 'RawAccel.csv'; 
gyroFile = 'RawGyro.csv'; 
GPSFile = 'OnboardGPS.csv'; 

accelTbl = readtable(accelFile); 
accelTimeStamp = accelTbl.Timpstemp;
accel_x       = accelTbl.x;
accel_y       = accelTbl.y;
accel_z       = accelTbl.z;

accelData = [accelTimeStamp, accel_x, accel_y, accel_z]; 

gyroTbl = readtable(gyroFile); 
gyroData = [gyroTbl.Timpstemp, gyroTbl.x, gyroTbl.y, gyroTbl.z]; 


gpsTbl = readtable(GPSFile); 
GPSData = [gpsTbl.Timpstemp,...
           gpsTbl.lat, gpsTbl.lon, gpsTbl.alt, ...
           gpsTbl.vel_n_m_s, gpsTbl.vel_e_m_s, gpsTbl.vel_d_m_s, ...
           gpsTbl.eph_m, gpsTbl.epv_m, ...
           gpsTbl.s_variance_m_s]; 


end