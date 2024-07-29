%% REACTION WHEELS 

% rw.a = 1/sqrt(3);
% rw.A = [-rw.a, rw.a, rw.a, -rw.a; -rw.a, -rw.a, rw.a, rw.a; rw.a, rw.a, rw.a, rw.a];
% rw.b = sqrt(3)/4;
% rw.Ainv = [-rw.b, -rw.b, rw.b; rw.b, -rw.b, rw.b; rw.b, rw.b, rw.b; -rw.b, rw.b, rw.b];

rw.A = [1, 0, 0, 1/sqrt(3); 0, 1, 0, 1/sqrt(3); 0, 0, 1, 1/sqrt(3)];
rw.Ainv = [5/6, -1/6, - 1/6; -1/6, 5/6, -1/6; -1/6, -1/6, 5/6; 1/(2*sqrt(3)), 1/(2*sqrt(3)), 1/(2*sqrt(3))];

% Saturation
rw.hMax = 1;
rw.dhMax = 0.1;
%% MAGNETORQUER

magnetorquer.m = 30; % A*m^2

%% INERTIA WHEEL

