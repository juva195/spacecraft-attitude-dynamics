%% EARTH ORBIT

mjd2000 = date2mjd2000(date);
[earth.kep, ~] = uplanet(mjd2000, 3); % Initial state of Earth's orbit in keplerian elements
[earth.r0, earth.v0] = orbitalToCar(earth.kep(1), earth.kep(2), earth.kep(3),... % Initial state of Earth's orbit in state vector
    earth.kep(4), earth.kep(5), earth.kep(6), astroConstants.muSun);

RAAN = atan(-earth.r0(1)/earth.r0(2));

earth.r0 = rotx(23.4)*earth.r0;
earth.v0 = rotx(23.4)*earth.v0;

%% SATELLITE ORBIT

orbit.a = ((86400/(15*2*pi))^2*astroConstants.muEarth)^(1/3); % Orbit semi-mayor axis [km] / 15 orbits per day
orbit.i = acos(-(orbit.a/12352)^(7/2)); % Orbit inclination [rad] / Sun-synchronous orbit
orbit.e = 0;
orbit.RAAN = RAAN;
orbit.omega = 0;
orbit.TrueAnomaly = 0;
orbit.T = 2 * pi * sqrt(orbit.a^3 / astroConstants.muEarth);

[orbit.r0, orbit.v0] = orbitalToCar(orbit.a,orbit.e,orbit.i,orbit.RAAN,orbit.omega,orbit.TrueAnomaly,astroConstants.muEarth);

track = [orbit.a,0,orbit.i,orbit.RAAN,0,0,2*pi];

%Ground track
default = [];
[~,~,~,~] = groundTrack(track,1,default,default,default,default,default);

%Orbit representation
Options = ["Ocompletion", "yes","setViasual","Elevation"];

orbitDraw([orbit.a,0,orbit.i,orbit.RAAN,0,0,2*pi],Options)

