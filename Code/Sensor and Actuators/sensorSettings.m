%% SUN SENSORS
sunSensor.dir = roty(90);
sunSensor.sampleRate = 2; % [Hz]
sunSensor.resolution = 0.2*pi/180; % [rad]
sunSensor.sd = sunSensor.resolution/(4); % [rad]

%% MAGNETOMETER
magSensor.sampleRate = 10; % [Hz]
magSensor.resolution = 8; % [nT]
magSensor.noisePower = 16^2; % @1Hz [nT^2 rms^2/Hz]

%% Gyroscope
% ASTRIX 200 Gyroscope

gyroscope.arw = 0.0001; % [°/(h)^0.5] 1sigma, found on the datasheet
% This is for Begging of Life value write in the report!!!
gyroscope.biasStability = 0.0005; % [°/h] 3sigma, found on the datasheet
gyroscope.f = 100; % [Hz] Samplinmg frequency, found on the datasheet

gyroscope.arw = gyroscope.arw/sqrt(3600); % [°/(s)^0.5]
gyroscope.biasStability = gyroscope.biasStability/(3*3600); % [°/s] 
gyroscope.Ts = 1/gyroscope.f; %[s] Sampling Time
gyroscope.stdDeviation = gyroscope.arw/gyroscope.Ts;