close all
clear all
clc

%% GENERAL SETTINGS AND CONSTANTS

sadSettings;

%% SENSORS SETTINGS

sensorSettings;

%% ACTUATORS

actuatorSettings;

%% INITIAL CONDITIONS

deploymentVelocity = 2; % Satellite deployment velocity [m/s]
deploymentMissalignment = 0.03; % Misalignment of deployment force [m]
satelliteMass = 133.2; % Satellite mass [kg]

h0 = [1;1;1] * deploymentVelocity * satelliteMass * deploymentMissalignment / sqrt(3);
w0 = h0./sat.I';
euler0 = [0; 0; 0]; % Initial orientation in 312 [rad]
date = [2023, 1, 21, 13, 52, 0]; % Mission day [year, month, day, hour, minute, second]

%% ORBIT DETERMINATION

OrbitDetermination

%% SIMULATION

movmeanWindow = 1000 * 1/(str2double(sim_options.FixedStep));

simout = sim('eulerEq_nonLin.slx', sim_options);

%% SIM OUTPUT

w = squeeze(simout.w);
dw = squeeze(simout.dw);
t = simout.t;
A = simout.A;
Adet = simout.Adet;
wdet = squeeze(simout.wdet);
hr = squeeze(simout.hr);
dhr = squeeze(simout.dhr);
disturbanceTorque = squeeze(simout.disturbanceTorque);
controlTorque = squeeze(simout.controlTorque);
disturbanceMomentum = squeeze(simout.disturbanceMomentum);
magneticTorque = squeeze(simout.magneticTorque);
sunVisible = squeeze(simout.sunVisible);
degPointingError = squeeze(simout.degPointingError);
state = squeeze(simout.state);
stateIntermediate = squeeze(simout.stateIntermediate);
derivative = squeeze(simout.derivative);
state2 = squeeze(simout.state2);
DCMTarget = simout.DCMTarget;
angle = simout.angle;

%% USEFUL VARIABLES

% Reaction wheel torque
dhrBody = rw.A * dhr;
dhrBody(1,:) = movmean(dhrBody(1,:), movmeanWindow);
dhrBody(2,:) = movmean(dhrBody(2,:), movmeanWindow);
dhrBody(3,:) = movmean(dhrBody(3,:), movmeanWindow);

% Reaction wheel momentum
hrBody = rw.A * hr;

% % Control state change
stateIdx = find(state>=1,1);
stateTime = t(stateIdx);

stateIdx2 = find(state2>=1,1);
stateTime2 = t(stateIdx2);


%% PLOTS


% Angular velocity and control state
f1 = figure();
plot(t,wdet)
xline(stateTime,'--','LineWidth',2)
xline(stateTime2,'--','LineWidth',2)
ylabel("$\omega$ [$\frac{rad}{s}$]",'Interpreter','latex')
yyaxis right
plot(t,movmean(squeeze(pagenorm(A-Adet)),100))
legend("$\omega_x$","$\omega_y$","$\omega_z$","slew state change","inertial pointing state change","DCM error")
title("$\omega$ measured",'Interpreter','latex')
ylabel("$movmean(||A - A_{measured}||)$",'Interpreter','latex')
xlabel("Time [$s$]",'Interpreter','latex')


f2 = figure();
plot(t,hr)
yline(-1,'LineWidth',2)
yline(1,'LineWidth',2)
legend("$h_r$ $1^{st}$ wheel","$h_r$ $2^{nd}$ wheel","$h_r$ $3^{rd}$ wheel","$h_r$ $4^{th}$ wheel")
title("Reaction Wheels Angular Momentum",'Interpreter','latex')
ylabel("$h_r$ [$Nms$]",'Interpreter','latex')
xlabel("Time [$s$]",'Interpreter','latex')

f3 = figure();
plot(t, degPointingError)
xline(stateTime,'--','LineWidth',2)
xline(stateTime2,'--','LineWidth',2)


% figure
% plot(t,i)
% hold on
% plot(t,state)
% legend("$\dot{\omega}_x$","$\dot{\omega}_y$","$\\dot{\omega}_z$","controller state")
% title("$\dot{\omega}$ measured",'Interpreter','latex')
% ylabel("$\dot{\omega}$ [$\frac{rad}{s^2}$]",'Interpreter','latex')
% xlabel("Time [$s$]",'Interpreter','latex')
% 
% % 
% tp = theaterPlot('XLimit',[-2 2],'YLimit',[-2 2],'ZLimit',[-2 2]);
% op = orientationPlotter(tp,'DisplayName','Fused Data',...
%     'LocalAxesLength',2);
% for i=1:length(A)
%     plotOrientation(op, cat(3,Adet(:,:,i),controldcmTarget))
%     drawnow
% end
% % 
% figure
%%
tp = theaterPlot('XLimit',[-2 2],'YLimit',[-2 2],'ZLimit',[-2 2]);
op = orientationPlotter(tp,'DisplayName','Fused Data',...
    'LocalAxesLength',2);
op2 = orientationPlotter(tp,'DisplayName','Fused Data',...
    'LocalAxesLength',2);
for i=1:length(DCMTarget)
    plotOrientation(op, DCMTarget(:,:,stateIdx2+i*10))
    drawnow
    hold on
    plotOrientation(op2, Adet(:,:,stateIdx2+i*10))
    title(stateIdx2+i*10)
    drawnow
end
