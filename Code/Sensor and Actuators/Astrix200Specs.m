% ASTRIX 200 Gyroscope

ARW = 0.0001; % deg/sqrt(hour) 1sigma

% This is for Begging of Life value write in the report!!!

BiasStability = 0.0005/3; % deg/hour 3sigma

ARW = ARW/sqrt(3600); % deg/sqrt(s) 1sigma
BiasStability = BiasStability/3600; %deg/s 3sigma

f = 100; %Hz
Ts = 1/f; %s

StdDevARW = ARW/Ts;



