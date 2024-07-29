function dy = ode_2bp( ~, y, mu, J, rE)
% ode_2bp ODE system for the two-body problem (Keplerian motion) both
% perturbated and not
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu, J, rE)
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% J[1] Second Zonal Harmonic [-]
% rE[1] Earth's radius [L]
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

if exist("y","var") && exist("mu","var")

    if exist("J","var") && exist("rE","var")
        % Set the derivatives of the state and add the perturbances
        aj2 = (3*J*mu*rE^2/(2*rnorm^4))*[(5*y(3)^2/(rnorm^2)-1)/rnorm (5*y(3)^2/(rnorm^2)-1)/rnorm (5*y(3)^2/(rnorm^2)-3)/rnorm]';
        dy = [ v; (-mu/rnorm^3)*r+aj2.*r ];
    else
        % Set the derivatives of the state withou peturbances
        dy = [ v; (-mu/rnorm^3)*r ];
    end
else
    error('y and/or mu are missing')
end


end