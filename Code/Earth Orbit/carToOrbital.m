function [a,eNorm,i,Omega,omega,theta,e] = carToOrbital(r,v,mu)
%
%% Conversion from Cartesian coordinates to Keplerian elements.           
% Angles are given in radians.
%
% PROTOTYPE:
% [a,e,i,Omega,omega,theta] = carToOrbital(r,v,mu)
%
% If default values are desired insert an empty vector ([]) as the corresponding value
%
% INPUT:
% r                 [3x1]           Position vector                                                                     [km]
% v                 [3x1]           Velocity vector                                                                     [km/s]
% mu                [1x1]           Gravitational parameter of the planet (for default mu = 3.98600433e+05)             [km^3/s^2]
%
% OUTPUT:
% a                 [1x1]           Semi-major axis                 [km]
% e                 [1x1]           Eccentricity                    [-]
% i                 [1x1]           Inclination                     [rad]
% Omega             [1x1]           RAAN                            [rad]
% omega             [1x1]           Pericentre anomaly              [rad]
% theta             [1x1]           True anomaly                    [rad]


%% VALUE CHECK
if nargin < 2
    error("Error: not enough variables")
end

if length(r) ~= 3 || length(v) ~= 3
    error("Please provide a 3x1 vector for both radius and velocity");
end

if isempty(mu)
    mu = astroConstants(13);
end

kDir = [0; 0; 1];
jDir = [0; 1; 0];
iDir = [1; 0; 0];

%% VECTOR NORMALIZATION

rNorm = norm(r);
vNorm = norm(v);

%% SEMI - MAJOR AXIS

a = 1/(2/rNorm - vNorm^2/mu);

%% ANGULAR MOMENTUM

h = cross(r,v);
hNorm = norm(h);

%% ECCENTRICITY

e = 1/mu * (cross(v,h) - mu*r/rNorm);
eNorm = norm(e);

%% INCLINATION

i = acos(dot(h,kDir)/hNorm);

%% NODAL PLANE DIRECTION

N = cross(kDir,h);
NNorm = norm(N);

%% RIGHT ASCENSION OF THE ASCENDING NODE - LONGITUDE OF ASCENDING NODE

if NNorm == 0
    Omega = 0;
else
    if dot(N,jDir) >= 0
        Omega = acos( dot(N,iDir) / NNorm);
    else 
        Omega = 2*pi - acos( dot(N,iDir) / NNorm);
    end
end

%% ARGUMENT OF PERIAPSIS

if NNorm == 0 || eNorm == 0
    omega = 0;
else
    if dot(e,kDir) >= 0 
        omega = acos(dot(N,e) / (NNorm * eNorm));
    else 
        omega = 2*pi - acos(dot(N,e) / (NNorm * eNorm));
    end
end

%% MEAN ANOMALY

if dot(v,r) >= 0
    theta = acos(dot(e,r) / (eNorm * rNorm));
else 
    theta = 2*pi - acos(dot(e,r) / (eNorm * rNorm));
end


end