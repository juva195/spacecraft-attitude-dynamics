%close all
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

simout = sim('eulerEq.slx', sim_options);

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
ggTorque = squeeze(simout.GG);
srpTorque = squeeze(simout.SRP);
mtTorque = squeeze(simout.MT);
adTorque = squeeze(simout.AD);
dipole = squeeze(simout.dipole);

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
sgtitle("Angular rates during de-tumbling")
subplot(3,1,1)
plot(t(1:stateIdx)/orbit.T,wdet(1,1:stateIdx))
xline(stateTime/orbit.T,'-','LineWidth',2)
xline(stateTime2/orbit.T,'--','LineWidth',2)
ylabel("$\omega \; \left[ 1/s\right]$",'Interpreter','latex')
legend("$\omega_x $","De-tumbling to non-linear control","Non-linear to linear control",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xlim([0, 2.1])


subplot(3,1,2)
plot(t(1:stateIdx)/orbit.T,wdet(2,1:stateIdx),'Color',[0.8500 0.3250 0.0980])
xline(stateTime/orbit.T,'-','LineWidth',2)
xline(stateTime2/orbit.T,'--','LineWidth',2)
ylabel("$\omega \; \left[ 1/s\right]$",'Interpreter','latex')
legend("$\omega_y $","De-tumbling to non-linear control","Non-linear to linear control",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xlim([0, 2.1])


subplot(3,1,3)
plot(t(1:stateIdx)/orbit.T,wdet(3,1:stateIdx),'Color',[0.9290 0.6940 0.1250])
xline(stateTime/orbit.T,'-','LineWidth',2)
xline(stateTime2/orbit.T,'--','LineWidth',2)
ylabel("$\omega \; \left[ 1/s\right]$",'Interpreter','latex')
legend("$\omega_z $","De-tumbling to non-linear control","Non-linear to linear control",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xlim([0, 2.1])



% Total angular momentum
f2 = figure();
plot(t(1:stateIdx)/orbit.T, sqrt((sat.Ixx * w(1,1:stateIdx)).^2 + (sat.Iyy * w(2,1:stateIdx)).^2 + (sat.Izz * w(3,1:stateIdx)).^2))
ylabel("Angular momentum $\left[ \frac{kg\;m^2}{s}\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xlim([0, 2.1])
title("Total angular momentum during de-tumbling",'Interpreter','latex')

% Control torque
f3 = figure();
plot(t(1:stateIdx)/orbit.T, vecnorm(controlTorque(1:stateIdx),2,1))
hold on
plot(t(1:stateIdx)/orbit.T, movmean(vecnorm(controlTorque(1:stateIdx),2,1),1000))
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
title("Control torque during de-tumbling",'Interpreter','latex')
legend("Norm of control torque","Average norm of control torque")
xlim([0, 2.1])

% Magnetic dipole
f4 = figure();
d1 = movmean(abs(dipole(1,:)),100);
d2 = movmean(abs(dipole(2,:)),100);
d3 = movmean(abs(dipole(3,:)),100);
plot(t(1:stateIdx)/orbit.T, d1(1:stateIdx))
hold on
plot(t(1:stateIdx)/orbit.T, d2(1:stateIdx))
plot(t(1:stateIdx)/orbit.T, d3(1:stateIdx))
yline(magnetorquer.m)
yline(-magnetorquer.m)
ylabel("Magnetic diple $\left[Am^2\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
title("Magnetorquers magnetic dipole during de-tumbling",'Interpreter','latex')
legend("$m_x$","$m_y$","$m_z$",'Interpreter','latex')
xlim([0, 2.1])
ylim([0, 30])

% Pointing error
f5 = figure();
plot(t(stateIdx:end)/orbit.T, degPointingError(stateIdx:end))
xline(stateTime2/orbit.T,'--','LineWidth',2)
ylabel("Error $\left[deg\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
title("Pointing error during slew and earth pointing",'Interpreter','latex')
legend("Pointing error", "Non-linear to linear control")
xlim([2.066,2.43])

% Wheels angular momentum
f6 = figure();
plot(t(stateIdx:end)/orbit.T,hr(:,stateIdx:end))
yline(-1,'-','LineWidth',2)
yline(1,'-','LineWidth',2)
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("$h_r$ $1^{st}$ wheel","$h_r$ $2^{nd}$ wheel","$h_r$ $3^{rd}$ wheel","$h_r$ $4^{th}$ wheel","Reaction wheels momentum saturation","","Non-linear to linear control")
title("Reaction Wheels Angular Momentum during slew and earth pointing",'Interpreter','latex')
ylabel("Angular momentum [$Nms$]",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xlim([2.066,2.43])
ylim([-1.1, 1.1])



% Control torque
f7 = figure();
plot(t(stateIdx:end)/orbit.T, vecnorm(controlTorque(:,stateIdx:end),2,1))
hold on
plot(t(stateIdx:end)/orbit.T, movmean(vecnorm(controlTorque(:,stateIdx:end),2,1),1000))
xline(stateTime2/orbit.T,'--','LineWidth',2)
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
title("Control torque during slew and earth pointing",'Interpreter','latex')
legend("Norm of control torque","Average norm of control torque")
xlim([2.066,2.43])

% Attitude error
f8 = figure();
plot(t/orbit.T, squeeze(pagenorm(A-Adet)))
plot(t/orbit.T,movmean(squeeze(pagenorm(A-Adet)),1000))
xline(stateTime/orbit.T,'-','LineWidth',2)
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("Attitude error", "Mean attitude error","De-tumbling to non-linear control","Non-linear to linear control")
title("Attitude error",'Interpreter','latex')
ylabel("Error $(||A - A_{measured}||)$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')

% Disturbance torques
f9 = figure();
sgtitle("Disturbance torques")
subplot(4,1,1)
plot(t(stateIdx:end)/orbit.T, movmean(vecnorm(ggTorque(stateIdx:end,:),2,2),1000))
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("Gravity gradient","Non-linear to linear control",'Location','best')

subplot(4,1,2)
plot(t(stateIdx:end)/orbit.T, movmean(vecnorm(srpTorque(:,stateIdx:end),2,1),1000))
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("Solar radiation pressure","Non-linear to linear control",'Location','best')

subplot(4,1,3)
plot(t(stateIdx:end)/orbit.T, movmean(vecnorm(mtTorque(stateIdx:end,:),2,2),1000))
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("Magnetic field","Non-linear to linear control",'Location','best')

subplot(4,1,4)
plot(t(stateIdx:end)/orbit.T, movmean(vecnorm(adTorque(:,stateIdx:end),2,1),1000))
ylabel("Torque $\left[ N\,m\right]$",'Interpreter','latex')
xlabel("Time [Orbits]",'Interpreter','latex')
xline(stateTime2/orbit.T,'--','LineWidth',2)
legend("Atmospheric drag","Non-linear to linear control",'Location','best')



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
