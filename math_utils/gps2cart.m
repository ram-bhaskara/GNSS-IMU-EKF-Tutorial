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
