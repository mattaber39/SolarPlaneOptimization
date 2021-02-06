function [P_thrust] = thrust_power_required(aero,geom,W,S,theta)
% Find thrust power required - code verified by hand calcs.
% Accounts for losses in propulsive efficiency.

% Extract variables
[~,~,~, rho] = atmosisa(aero.alt);
mu = 1.9*10^-5;          %dynamic viscosity @ 40 C

% Find thrust lapse (alpha) due to altitude
[~,~,~, rho0] = atmosisa(0);
alpha = rho/rho0;
    
%% Use drag polar to find C_L and C_D

% P_given = 100;
[C_L, C_D, ~] = Drag_Polar(aero,geom,S,W);   %C_L and C_D will already be minimum power required values
% v = sqrt(2*W*cos(theta)/(C_L*S*rho));
% 
% Drag   = 0.5*rho*S*C_D*v^2;
% theta_fall = asin(-C_D*cos(theta)/(C_L));
% theta_fallView = theta_fall*180/pi;
% 
% theta_climb = asin(P_given/v/W-C_D*cos(theta)/(C_L));
% theta_climbView = theta_climb*180/pi;
% 
% Thrust_fall = Drag + W*sin(theta_fall);
% P_thrust_fall = Thrust_fall*v/aero.eta_prop;
% 
% Thrust_climb = Drag+W*sin(theta_climb);
% P_thrust_fall = Thrust_fall*v/aero.eta_prop;

% P_thrust_climb = 
P_thrust = (2/(rho*S))^.5*(W/C_L)^1.5*C_D/(alpha*aero.eta_prop);
%% calculate Re number
% chord = sqrt(S/aero.AR)/2;         %SMC
% Re = (rho*v*chord)/mu;

%find C_L that A/C is flying at in SLF
% C_L = W/(.5*rho*v^2*S);

%% Find thrust power required
% Derived from P = Drag*v, where 
% v is found using Lift = Weight = 0.5*rho*S*Cl*v^2, and 
% Drag = 0.5*rho*S*Cd*v^2. 
end
