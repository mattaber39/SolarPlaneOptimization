function [P_thrust] = thrust_power_required(aero,geom,W,S)
% Find thrust power required - code verified by hand calcs.
% Accounts for losses in propulsive efficiency.

% Extract variables
[~,~,~, rho] = atmosisa(aero.alt);

% Find thrust lapse (alpha) due to altitude
[~,~,~, rho0] = atmosisa(0);
alpha = rho/rho0;

% Use drag polar to find C_L and C_D
[C_L, C_D] = Drag_Polar(aero,geom,S,W);   %C_L and C_D will already be minimum power required values

% Find thrust power required. Derived from P = Drag*v, where 
% v is found using Lift = Weight = 0.5*rho*S*Cl*v^2, and 
% Drag = 0.5*rho*S*Cd*v^2. 
P_thrust = (2/(rho*S))^.5*(W/C_L)^1.5*C_D/(alpha*aero.eta_prop);
end

