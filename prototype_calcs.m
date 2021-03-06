function [n_s] = prototype_calcs(span, AR)
% README:
% This is a program to quickly calculate some important features of the 
% first prototype, such as # of solar panels, voltage, current, power, and 
% energy produced by the solar panels. There is no optimization functionality 
% because the focus is on testing the electronics on a platform easy to build. 

% Assumptions:
% The panels on each wing are wired in series and then each wing is wired in
% parellel with each other. 
% The worst specs were used (Ultra High Performance) for the solar panels.

% User input parameters
% span = 5;       % Wingspan [m]
% AR = 20;        % Aspect Ratio

% Calculated/User Input Parameters
chrd = span/AR;     % Chord [m]
wl = span/3;        % Length of one wing
cs_w = wl;          % Control Surface width on each side [m]
cs_d = chrd/4;      % Control Surface depth into the wing [m]
max_cover = 0.75;   % Amount of airfoil able to be covered by solar panels

% Technological parameters
sp_l = 0.125;   % Solar panel length and width [m]
sp_A = sp_l^2;  % Solar panel area [m^2]
Vmpp = 0.61;    % Voltage produced per cell at max power point [V]
Voc  = 0.71;    % Max voltage each cell is capable of [V]
Impp = 5.8;     % Current produced per cell at max power point [A]
Isc  = 6.1;     % Max current each cell is capable of [A]

% Calculations
rows = floor(chrd*max_cover/sp_l);               % # of rows for non-control-surface sections.
cs_rows = floor((chrd*max_cover - cs_d)/sp_l);   % # of rows for control surfaces sections.
clmns = floor((wl - cs_w)/sp_l);                 % # of columns for non-control-surface sections.
cs_clmns = floor(cs_w/sp_l);                     % # of columns for control surface sections.
covered = (100*rows*sp_l/chrd);

N = (rows*clmns+cs_rows*cs_clmns);   % # of panels per wing.

Vmpp_tot = N*Vmpp;  % Voltage each wing produces at max power point [V].
Voc_tot  = N*Voc;   % Max voltage each wing is capable of [V].

Pmpp = Vmpp_tot*Impp*2;   % Total power output per at max power point [kW].
% E_tot = 
 
n_s = N*2;
% fprintf("\nUser Inputs\n");
% fprintf("Wingspan: %d m \n", span);
% fprintf("Aspect Ratio: %d \n", AR);
% fprintf("Amount of Airfoil Covered: %d percent \n\n", covered);
% disp("Main info");
% fprintf("# of Panels per Wing: %d\n", N);
% fprintf("Total Voltage at max power point: %.1f V \n", Vmpp_tot);
% fprintf("Total Current at max power point: %.1f A \n", Impp*2);
% fprintf("Total Power   at max power point: %.1f W \n\n", Pmpp);
% disp("Additional Info");
% fprintf("Total Peak Voltage: %.1f V \n", Voc_tot);
% fprintf("Total Peak Current: %.1f A \n", Isc*2);


