function [r,v] = orbitalToCar(a,e,i,Omega,omega,theta,mu)
%
% Conversion from Keplerian elements to Cartesian coordinates.
% Angles must be given in radians.
%
% PROTOTYPE:
% [r, v] = orbitalToCar(a,e,i,Omega,omega,theta)
%
% If default values are desired insert an empty vector ([]) as the corresponding value
%
% INPUT:
% a         [1x1]         Semi-major axis                                                                     [km]
% e         [1x1]-[3x1]   Eccentricity                                                                        [-]
% i         [1x1]         Inclination                                                                         [rad]
% Omega     [1x1]         RAAN                                                                                [rad]
% omega     [1x1]         Pericentre anomaly                                                                  [rad]
% theta     [1x1]         True anomaly                                                                        [rad]
% mu        [1x1]         Gravitational parameter of the planet (for default mu = 3.98600433e+05)             [km^3/s^2]
%
%
% OUTPUT:
% r         [3x1]   Position vector [km]
% v         [3x1]   Velocity vector [km/s]

%% VALUE CHECK
if nargin < 6
    error("Please insert a valid amount of variables");
end

if length(e) == 3 
    e = norm(e);
elseif length(e) ~= 1
    error("Please check e value, it is neither a 3 element vector nor a scalar!");
end

if isempty(mu)
    mu = astroConstants(13);
end

%% R,V PERIFOCAL SYSTEM
if e >= 0 && e < 1 % Elliptic
    p = a*(1-e^2); % Semilat rett
elseif e == 1 % Parabolic
    p = 2 * a;
else
    p = a*(e^2 - 1);    
end

rNorm = p/(1+e*cos(theta)); % Radius
rPF = rNorm * [cos(theta); sin(theta); 0];
vPF = sqrt(mu/p) * [-sin(theta); (e+cos(theta)); 0];
%% FROM PERIFOCAL TO GEOCENTRIC

% Rotation of angle Omega along earth's K axis
R1 = [cos(Omega),sin(Omega),0; -sin(Omega),cos(Omega),0; 0,0,1];

% Rotation of angle i along our orbit's I axis
R2 = [1,0,0; 0,cos(i),sin(i); 0,-sin(i),cos(i)];

% Rotation of angle omega along our orbit's K axis
R3 = [cos(omega),sin(omega),0; -sin(omega),cos(omega),0; 0,0,1];

% Rotation
T = R1'*R2'*R3';

% From Perifocal to Geocentric Equatorial
r = T * rPF;
v = T * vPF;

end

