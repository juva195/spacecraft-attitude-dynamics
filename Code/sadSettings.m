%% PATH 

if ismac()
    addpath(genpath("./"));
else
    addpath(genpath(".\"));
end
%% CONSTANTS 

% Astronomical Constants
astroConstants.solRadP = 1358; %Solar radiation intensity [W/m^2]
astroConstants.aEarth = 6378.137; % Semi major axis of earth (a) [km]
astroConstants.bEarth = 6356.75231414; % Semi minor axis of earth (b) [km]
astroConstants.wEarth = 0.000072921; % Angular velocity of earth [rad/s]
astroConstants.muEarth = 3.986004418e5; % Standard gravitational parameter of Earth [km^3/s^2]
astroConstants.muSun = 1.32712440018e11; % Standard gravitational parameter of the Sun [km^3/s^2]

% Physics constants
physicsConstants.Cd = 2; % Coeficient of drag of a flat plane
physicsConstants.c = 299792458; % Speed of light [m/s]
airDensity.baseAltitude = [300, 350, 400, 450, 500, 600, 700, 800, 900, 1000]; % [km]
airDensity.nominalDensity = [2.418e-11, 9.158e-12, 3.725e-12, 1.585e-12, 6.967e-13, 1.454e-13, 3.614e-14, 1.17e-14, 5.245e-15, 3.019e-15]; % [kg/m^3]
airDensity.scaleHeight = [53.628, 53.298, 58.515, 60.828, 63.822, 71.835, 88.667, 124.64, 181.05, 268]; % [km]

%%  SATELLITE GEOMETRIC AND PHYSICAL PROPERTIES

sat.NB = [[1;0;0],[0;1;0],[-1;0;0],[0;-1;0],[0;0;1],[0;0;-1],[0;1;0],[0;-1;0]]; % Face normals
sat.ps = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.1]; % Specular coeficient
sat.pd = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; % Difusive coeficient
sat.A = [0.48, 0.48, 0.48, 0.48, 0.36, 0.36, 0.72, 0.72]; % Face areas [m^2]
sat.rF = [[30;0;0],[0;30;0],[-30;0;0],[0;-30;0],[0;0;37.2973],[0;0;-42.7027],[0;0;97.2973],[0;0;97.2973]]*10^-2; % Face center of pressure [m]
sat.Ixx = 14842702.702703000 * 1e-6;
sat.Iyy = 14734732.702703000 * 1e-6;
sat.Izz = 7884030.000000000 * 1e-6;
sat.I = [sat.Ixx, sat.Iyy, sat.Izz]; % Moment of inertia [kg*m^2]
sat.m = [0.1332;0.1332;0.1332]; % Magnetic dipole [A/m^2]

%% MAGNETIC FIELD 

[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
magneticConstants.igrfSg = [gn, gm, gvali, gsvi];
magneticConstants.igrfSh = [hn, hm, hvali, hsvi];

%% CONTROL ALGORITHM SETTINGS

control.wMinDetumbling = 0.01;
control.pointingErrorMaxLinear = 3;
control.nominalHw = [0, 0, 0];
control.dcmTarget = [-0.7574, 0.5143, -0.4024; -0.3758, -0.8472, -0.3756; -0.5341, -0.1333, 0.8349];
control.qTarget = dcm2quat(control.dcmTarget);
%control.qTarget = [1,1,1,1]/2;
control.wTarget = [0; 2*pi/5760; 0];
%% STATE SPACE REPRESENTATION

state.A = [ 0, 0, (sat.Izz-sat.Iyy)/sat.Ixx * control.wTarget(2) - control.nominalHw(2), 0, 0, 0, 0;...
            0, 0, 0, 0, 0, 0, 0;...
            (sat.Iyy-sat.Ixx)/sat.Izz * control.wTarget(2) + control.nominalHw(2), 0, 0, 0, 0, 0, 0;...
            -1/2 * control.qTarget(2), -1/2 * control.qTarget(3), -1/2 * control.qTarget(4), 0, 0, 0, 0;...
            1/2 * control.qTarget(1), -1/2 * control.qTarget(4), 1/2 * control.qTarget(3), 0, 0, 0, 0;...
            1/2 * control.qTarget(4), 1/2 * control.qTarget(1), -1/2 * control.qTarget(2), 0, 0, 0, 0;...
            -1/2 * control.qTarget(3), 1/2 * control.qTarget(2), 1/2 * control.qTarget(1), 0, 0, 0, 0];
state.B = [ 1/sat.Ixx, 0, 0;...
            0, 1/sat.Iyy, 0;...
            0, 0, 1/sat.Izz;...
            0, 0, 0;...
            0, 0, 0;...
            0, 0, 0;...
            0, 0, 0];

state.R = diag(1 * ones(1,3));
state.Q = diag([100, 100, 100, 1e-4, 1e-4, 1e-4, 1e-4]);
control.K = lqr(state.A, state.B, state.Q, state.R);

%% SIMULATION SETTINGS

sim_options.SolverType = 'Fixed-step';      % Set the solver type to Fixed-step
sim_options.Solver = 'ode4';                % Select ode4 as solver
sim_options.FixedStep = '0.01';              % Select timestep
sim_options.StartTime = '0';                % Start from 0 seconds [default]
sim_options.StopTime = '14000';                % End the simulation at

%% PLOT SETTINGS

set(groot, 'defaultTextInterpreter', 'latex')
set(groot,'defaultAxesXMinorGrid','on','defaultAxesXMinorGridMode','manual');
set(groot,'defaultAxesYMinorGrid','on','defaultAxesYMinorGridMode','manual');
set(groot, 'defaultLegendLocation', 'northeast')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter', 'latex')
set(groot, 'defaultAxesFontWeight', 'bold')
set(groot, 'defaultFigurePosition', [470, 360, 900, 530])
set(groot, 'defaultFigureColormap', turbo(256));
set(groot, 'defaultAxesFontName', 'Palatino Linotype', 'defaultTextFontName', 'Palatino Linotype');
set(groot, 'defaultSurfaceEdgeAlpha', 0.3);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultFigureColor', [1; 1; 1]);
set(groot, 'defaultAxesColor', 'none');
set(groot, 'defaultAxesFontSize', 20);
