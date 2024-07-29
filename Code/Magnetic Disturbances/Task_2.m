clear all

I = [0.07; 0.055; 0.025];
w0 = [0.45; 0.52; 0.55];
A0 = eye(3);
eulerA0 = [0; 0; 0];
q0 = [0; 0; 0; 1];
t0 = num2str(0);
tf = num2str(0.1);

%% Important stuff for air drag

astroConstants.c = 299792458; % Speed of light [m/s]
astroConstants.solRadP = 1358; % Solar radiation intensity [W/m^2]
astroConstants.aEarth = 6378137; % Semi major axis of earth (a) [m]
astroConstants.bEarth = 6356752.31414; % Semi minor axis of earth (b) [m]
astroConstants.wEarth = 0.000072921; % Angular velocity of earth [rad/s]

physicsConstants.Cd = 2; % Coeficient of drag of a flat plane

airDensity.baseAltitude = [300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]*10^3; % [m]
airDensity.nominalDensity = [2.418e-11, 9.158e-12, 3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.17e-14, 5.245e-15, 3.019e-15]; % [kg/m^3]
airDensity.scaleHeight = [53.628, 53.298, 58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268]*10^3; % [m]

% Satelite shape properties

sat.NB = [[1;0;0],[0;1;0],[-1;0;0],[0;-1;0],[0;0;1],[0;0;-1],[1;0;0],[-1;0;0],[1;0;0],[-1;0;0]]; % Face normals

sat.ps = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1, 0.1]; % Specular coeficient

sat.pd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; % Difusive coeficient

sat.A = [6, 6, 6, 6, 4, 4, 12, 12, 12, 12]*10^-2; % Face areas [m^2]

sat.rF = [[10;0;0],[0;10;0],[-10;0;0],[0;-10;0],[0;0;15],[0;0;-15],[0;0;45],[0;0;45],[0;0;-45],[0;0;-45]]*10^-2; % Face center of pressure [m]

%% Simulation
        
sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.1';              % Select a time step of -0.9 s
sim_options.StartTime = t0;                 % Start from 0 seconds [default]
sim_options.StopTime = tf;                  % End the simulation at t=tf

simout = sim('airDrag.slx', sim_options);

