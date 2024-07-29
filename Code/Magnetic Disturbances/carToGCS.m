function [rNorm, lat, lon] = carToGCS(r, T)
rNorm = norm(r);

omegaE = 15.04*pi/(3600*180);
thetaG0 = 0;
t0 = date2mjd2000([2023, 3, 20, 21, 24, 0]);

%% Define Latitude and Longitude

delta = asin(r(3)/rNorm);
alpha = atan2(r(2),r(1));

thetaG = omegaE*(T-t0) + thetaG0;

lon = rad2deg(wrapTo2Pi(alpha - thetaG + pi) - pi);         % Longitude
lat = rad2deg(delta);                                       % Latitude
end