function OM = quatOmega(omega)
    wx = omega(1); wy = omega(2); wz = omega(3);
    OM = [  0   -wx -wy -wz;
           wx    0   wz -wy;
           wy  -wz   0   wx;
           wz   wy -wx   0];
end