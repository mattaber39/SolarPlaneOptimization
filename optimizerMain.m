%% Setting things up and getting user input
clc; clear all; close all;
format compact
optimization = input("Enter 1 to run optimization\n");     
sensitivity = input("Enter 1 to run sensitivity study\n");

%% Initializing Power Train Info
% Assuming Samsung 36G 18650s and the solar panels we bought. 

PT.SunData    = 'SunData_060618';% String of text file name of sun data ('SunData_mmddyy')
PT.array_vert = 0;
PT.array_norm = 90;
PT.A          = 0.125^2;        % Area of each solar panel [m2]
PT.sigma      = 7;              % Weight per area of solar cells (measured from current array)   [N/m^2];
PT.W_s = 9.6;                   % Estimated total weight of solar array, assuming 150 panels.[N]
PT.eta_solar  = 1.0;            % Solar cell efficiency (from testing). The 
                                % efficiency has already been taken into
                                % account in the Solar_Array_Model
                                % function, when initializing I_nom.

PT.eta_mppt = .95;              % Efficiency of Max Power Point Tracker.
PT.eta_batt = .95;              % Efficiency of battery (will need further testing to validate).
PT.rho_batt = 270;              % Specific energy of battery (will need further testing to validate) [Wh/kg]
PT.batt_range = 0.8;            % Fraction of the batteries capacity that will be used (cannot have full discharge)
PT.e_b_raw = 46656;             % Raw energy of each battery [J].
PT.e_b = PT.e_b_raw*PT.batt_range; % Available Energy. 
PT.m_b = 0.048;                % Mass of each battery cell.
PT.w_b = PT.m_b*9.81;           % Weight of each battery cell.
% PT.eta_cc = .95;

%% Initializing Aero Info
% Possible airfoils: SA7024, MA409. Still using Luke's aero info. 

geom.t_c = 0.15;% Wing thickness/chord length

aero.alt = 100;                 %Altitude [m]

aero.C_DP = .0050;              %parasitic drag EXCLUDING contribution from wings
aero.h = .2;                    %vertical aspect ratio. Wing separation height/span

aero.k_dp = .005;               %k" value for airfoil. Used for lift-dependent viscous drag approximation 

aero.e = 1.3;                   %spanwise efficiency
aero.C_L_max = 1.1;             %max C_L
aero.eta_prop = 0.60;           %propulsive efficiency - will eventually be replaced with an engine deck.

%% Optimization Script
% Change the independant variables. 

% Independant Variables
n_b     = 50:2:80;
b       = 5:7;
AR      = 7:20;

% Finding sizes of independant variable arrays. 
n_size       = length(n_b);
b_size       = length(b);
AR_size      = length(AR);

% Initializing counters. 
n_count  = 1;
b_count  = 1;
AR_count = 1;

% Initializing storage for intermediate optimization steps. 
n_for_AR_loop = zeros(1, AR_size);
b_for_AR_loop = zeros(1, AR_size);
n_index_4b = zeros(1, b_size);
n_index_4AR = zeros(1, b_size);
b_index_4AR = zeros(1, b_size);
t_batt_AR     = zeros(1, AR_size);

% This part is uncommented and is not really in a good place to be read or
% understood.
if optimization == 1
for AR_count = 1:AR_size
    t_batt_b = zeros(1, b_size);
    for b_count = 1:b_size
        t_batt_n = zeros(1,n_size);
        for n_count = 1:n_size
            this_n_b = n_b(n_count);
            this_b = b(b_count);
            aero.AR = AR(AR_count);
            [t_batt, t_rech_night, t_rech_batt] = function2Optimize(this_n_b, this_b, aero.AR, aero, PT, geom);
            t_batt_n(n_count) = t_batt;  % Optimizing metric. 
        end
        [this_n_maxt, max_n_t_index] = max(t_batt_n);
        n_index_4b(b_count) = n_b(max_n_t_index);
        t_batt_b(b_count) = this_n_maxt;
    end
    [this_b_maxt, max_b_t_index] = max(t_batt_b);
    n_index_4AR(AR_count) = n_index_4b(max_b_t_index);
    b_index_4AR(AR_count) = b(max_b_t_index);
    t_batt_AR(AR_count) = this_b_maxt;
end

[this_AR_maxt, max_AR_t_index] = max(t_batt_AR);

n_b = n_index_4AR(max_AR_t_index);
b = b_index_4AR(max_AR_t_index);
aero.AR = AR(max_AR_t_index);

end
%% Non-Optimization/Output Step

% If optimization is not run, the inputs must be initialized.
if optimization ~= 1
n_b     = 60;           % Number of batteries
b       = 7;            % Span
aero.AR = 15;           % Aspect ratio
end

% Find outputs at the desired inputs. 
[t_batt, t_rech_night, t_rech_batt, t_day, t_night...
    n_s, M_total, M_e, v] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);

% Display outputs. 
fprintf("----------------------Inputs-----------------------\n")
fprintf("Number of batteries: %d\n", n_b)
fprintf("Wing span..........: %d\n", b)
fprintf("Aspect Ratio.......: %d\n", aero.AR)
fprintf("-----------------------Misc------------------------\n")
fprintf("Total mass..............: %.1f kg\n", M_total)
fprintf("Empty mass..............: %.1f kg\n", M_e)
fprintf("Battery mass............: %.1f kg\n", M_total-M_e)
fprintf("Battery/Total Mass Ratio: %.2f\n", (M_total-M_e)/M_total)
fprintf("Speed...................: %.1f m/s\n", v)
fprintf("Number of Solar Panels..: %d\n", n_s)
fprintf("-----------------------Outputs---------------------\n")
fprintf("Flight time with only battery: %.1f hrs\n", t_batt)
fprintf("Time to recharge for night...: %.1f hrs\n", t_rech_night)
fprintf("Time to fully recharge.......: %.1f hrs\n", t_rech_batt)
fprintf("-----------------------Margins---------------------\n")
fprintf("Extra flight time with only battery......: %.1f hrs\n", t_batt-t_night)
fprintf("Extra daylight after recharging for night: %.1f hrs\n", t_day-t_rech_night)
fprintf("Extra daylight after fully recharging....: %.1f hrs\n", t_day-t_rech_batt)


%% Plotting Sensitivities
% This part is uncommented and is not really in a good place to be read or
% understood. 

if sensitivity == 1
    
% Time of year:
this_SunData = PT.SunData;

SunData = {'SunData_030119', 'SunData_040119', 'SunData_060618', 'SunData_122119'};
dates   = [3 4 6 12];
SD_size = length(SunData);
SD_t_batt = zeros(1, SD_size);
SD_t_night = zeros(1, SD_size);
x_axis = zeros(1, SD_size);
for SD_count = 1:SD_size
    this_SD = SunData(SD_count);
    PT.SunData    = this_SD{:};
    [SD_t_batt(SD_count), ~, ~, ~, SD_t_night(SD_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
SD_t_extra = SD_t_batt - SD_t_night;
subplot(2, 3, 1);
plot(dates, SD_t_extra, 'x', 'Color', [0, 0.4, 0.7], 'MarkerSize', 10);
hold on;
plot(dates, x_axis, 'r')
plot(dates(3), SD_t_extra(3), 'og', 'MarkerSize', 10);
title("Time of Year Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Month of Year");
PT.SunData = 'SunData_060618';  % Return SunData to original value. 

% Altitude:
this_Altitude = 100;
aero.alt = this_Altitude;
[this_alt_t_batt, ~, ~, ~, this_t_night] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
this_alt_t_extra = this_alt_t_batt - this_t_night;
this_Altitude = this_Altitude*3.2;
Altitude = 50:10:200;
Alt_size = length(Altitude);
Alt_t_batt = zeros(1, Alt_size);
Alt_t_night = zeros(1, Alt_size);
for sc_count = 1:Alt_size
    aero.alt = Altitude(sc_count);
    [Alt_t_batt(sc_count), ~, ~, ~, Alt_t_night(sc_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
Alt_t_extra = Alt_t_batt - Alt_t_night;
Altitude = Altitude*3.2;
subplot(2, 3, 2)
plot(Altitude, Alt_t_extra);
hold on
plot(this_Altitude, this_alt_t_extra, 'og', 'MarkerSize', 10);
title("Altitude Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Altitude [ft]");
aero.alt = 100; % Return Altitude to original value. 

% Propeller Efficiency:
this_eta_prop = aero.eta_prop;
[this_ep_t_batt, ~, ~, ~, this_t_night] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
this_br_t_extra = this_ep_t_batt - this_t_night;
eta_prop = 0.1:0.05:1.0;
ep_size = length(eta_prop);
ep_t_batt = zeros(1, ep_size);
ep_t_night = zeros(1, ep_size);
x_axis = zeros(1, ep_size);
for sc_count = 1:ep_size
    aero.eta_prop = eta_prop(sc_count);
    [ep_t_batt(sc_count), ~, ~, ~, ep_t_night(sc_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
sc_t_extra = ep_t_batt - ep_t_night;
subplot(2, 3, 3)
plot(eta_prop, sc_t_extra, eta_prop, x_axis, 'r');
hold on
plot(this_eta_prop, this_br_t_extra, 'og', 'MarkerSize', 10);
title("Propeller Efficiency Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Efficiency");
legend("Sensitivity Study", "Lowest Allowable Flight Time", "Current Position")
newPosition = [0.4 0.4 0.2 0.2];
aero.eta_prop = 0.6; % Return to original value. 

% Battery Efficiency:
this_eta_batt = PT.eta_batt;
[this_batt_t_batt, ~, ~, ~, this_t_night] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
this_br_t_extra = this_batt_t_batt - this_t_night;
eta_batt = 0.6:0.05:1.0;
batt_size = length(eta_batt);
batt_t_batt = zeros(1, batt_size);
batt_t_night = zeros(1, batt_size);
x_axis = zeros(1, batt_size);
for Alt_count = 1:batt_size
    PT.eta_batt = eta_batt(Alt_count);
    [batt_t_batt(Alt_count), ~, ~, ~, batt_t_night(Alt_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
ep_t_extra = batt_t_batt - batt_t_night;
subplot(2, 3, 4)
plot(eta_batt, ep_t_extra, eta_batt, x_axis, 'r');
hold on
plot(this_eta_batt, this_br_t_extra, 'og', 'MarkerSize', 10);
title("Battery  Efficiency Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Efficiency");
PT.eta_batt = 0.95; % Return to original value. 

% Battery Discharge Depth:
this_batt_range = 0.8;
PT.e_b = PT.e_b_raw*PT.batt_range; % Available Energy. 
PT.batt_range = this_batt_range;
[this_br_t_batt, ~, ~, ~, this_t_night] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
this_br_t_extra = this_br_t_batt - this_t_night;
batt_range = 0.6:0.05:1.0;
br_size = length(batt_range);
br_t_batt = zeros(1, br_size);
br_t_night = zeros(1, br_size);
x_axis = zeros(1, br_size);
for br_count = 1:br_size
    PT.batt_range = batt_range(br_count);
    PT.e_b = PT.e_b_raw*PT.batt_range; % Available Energy. 
    [br_t_batt(br_count), ~, ~, ~, br_t_night(br_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
br_t_extra = br_t_batt - br_t_night;
subplot(2, 3, 5)
plot(batt_range, br_t_extra, batt_range, x_axis, 'r');
hold on
plot(this_batt_range, this_br_t_extra, 'og', 'MarkerSize', 10);
title("Battery  Discharge Depth Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Discharge Depth");
PT.batt_range = this_batt_range; % Return to original value. 
PT.e_b = PT.e_b_raw*PT.batt_range; % Available Energy. 

% Mass Sensitivity
this_mb = PT.m_b;
PT.w_b = PT.m_b*9.81;
[this_mb_t_batt, ~, ~, ~, this_mb_t_night, ~, this_M_total] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
this_t_extra = this_mb_t_batt - this_mb_t_night;
mb = 0.04:0.005:0.11;
mb_size = length(mb);
mb_t_batt = zeros(1, mb_size);
mb_t_night = zeros(1, mb_size);
mb_M_total = zeros(1, mb_size);
x_axis = zeros(1, mb_size);
for mb_count = 1:mb_size
    PT.m_b = mb(mb_count);
    PT.w_b = PT.m_b*9.81;           % Weight of each battery cell.
    [mb_t_batt(mb_count), ~, ~, ~, mb_t_night(mb_count), ~, mb_M_total(mb_count)] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
end
mb_t_extra = mb_t_batt - mb_t_night;
subplot(2, 3, 6)
plot(mb_M_total, mb_t_extra, mb_M_total, x_axis, 'r');
hold on
plot(this_M_total, this_t_extra, 'og', 'MarkerSize', 10);
title("Mass Sensitivity");
ylabel("Extra Battery Flight Time [hrs]"); xlabel("Total Mass [kg]");
PT.m_b = this_mb; % Return to original value. 
PT.w_b = PT.m_b*9.81;


% % Solar Cell Efficiency:
% eta_sc = 0.0:0.1:1.0;
% sc_size = length(eta_sc);
% sc_t_rech = zeros(1, sc_size);
% sc_t_day = zeros(1, sc_size);
% x_axis = zeros(1, sc_size);
% for sc_count = 1:sc_size
%     PT.eta_solar = eta_sc(sc_count);
%     [~, sc_t_rech(sc_count), ~, sc_t_day(sc_count), ~] = function2Optimize(n_b, b, aero.AR, aero, PT, geom);
% end
% sc_t_day
% sc_t_rech
% sc_t_extra = sc_t_day - sc_t_rech
% figure(5)
% plot(eta_sc, sc_t_rech, eta_sc, x_axis, 'r');
% title("Solar  Efficiency Sensitivity");
% ylabel("Extra Battery Flight Time [hrs]"); xlabel("Efficiency");
% PT.eta_solar = 1.0; % Return to original value. 

end
%% Function To Optimize
function [t_batt, t_rech_night, t_rech_batt, t_day, t_night ...
    n_s, M_total, M_e, v] = function2Optimize(n_b, b, AR, aero, PT, geom)
% Still messy. Will clean up soon. 

aero.AR = AR;
S = b^2/aero.AR;            % Area [m2].
E_b = n_b*PT.e_b;           % Total energy of battery pack [J].
W_b = n_b*PT.w_b;           % Find weight of batteries [N].
W_e = weight_model...       % Find empty weight using Lukes' weight model [N].
    (0,b,0,0,geom,aero,0);
W_t = W_e+W_b+PT.W_s;       % Find total estimated weight [N].
                            % PT.W_s = estimated weight of solar panels.
                            
% Find thrust required. 
P_flight = thrust_power_required(aero,geom,W_t,S);

% Find e_s, the energy each solar panel creates over the course of the day,
% in [W]. Find t_day, t_night (the length of day and night), in seconds. 
% Find the vector of power given by each solar panel throughout the day. 
[e_s, p_vec, t_day, t_night, ~, ts] = Solar_Array_Model(PT);

% Find the number of solar cells able to fit on the wing. 
n_s = prototype_calcs(b, aero.AR);    

% Find power vector given by the whole solar array throughout the day. 
P_vec = n_s*p_vec;

% Find residual of power required to fly after subtracting power input by
% the solar array. 
P_flight_res = P_flight - P_vec;

% Find the start of dusk. 
p_count = floor(length(P_vec)/2);
dusk_start_index = 0;
while P_flight_res(p_count) < 0
    dusk_start_index = p_count + 1;
    p_count = p_count + 1;
end

% Finding energy required to fly through dawn
Emech_dawn = 0; 
t_dawn = 0;     
p_count = 1;
while P_flight_res(p_count) > 0 && p_count < length(P_vec)
    Emech_dawn = Emech_dawn + P_flight_res(p_count)*ts;
    t_dawn = t_dawn + ts;
    p_count = p_count + 1;
end

% Finding energy required to fly through dusk
Emech_dusk = 0;
t_dusk = 0;
p_count = length(P_vec);
while P_flight_res(p_count) > 0 && p_count > 1
    Emech_dusk = Emech_dusk + P_flight_res(p_count)*ts;
    t_dusk = t_dusk + ts;
    p_count = p_count - 1;
end

% Finding energies
Emech_night = P_flight*t_night;           % Thrust energy need to fly during night.
Echem_dawn  = Emech_dawn/PT.eta_batt;     % Battery energy needed to fly during dawn.
Echem_dusk  = Emech_dusk/PT.eta_batt;     % Battery energy needed to fly during dusk.
Echem_night = Emech_night/PT.eta_batt;    % Battery energy needed to fly during night.

% Finding how long it takes to recharge for night. Assuming worst case
% scenario (starting at the very beginning of the day).
B_level = 0;            % Energy level of battery pack [J]
t_rech_night = 0;       % time to recharge for night only. 
p_count = t_dawn/ts+1;  % Start when solar array is giving excess power. 
while B_level < Echem_dusk + Echem_night + Echem_dawn && p_count <= dusk_start_index % While battery isn't charged for night and dusk hasn't arrived:
    B_level = B_level - P_flight_res(p_count)*ts;  % Add energy gained during 'ts' time to battery level. Since P_flight_res is negative, subtract it. 
    t_rech_night = t_rech_night + ts;  % Increment time. 
    if p_count == dusk_start_index
%         disp("Batteries did not recharge to last the next night cycle")
    end
    p_count = p_count + 1;      % Increment p-counter
end

% Finding out how long it takes to fully recharge.
t_rech_batt  = t_rech_night;    % Start where we left off. 
while B_level < E_b && p_count <= dusk_start_index % While battery isn't full and dusk hasn't arrived:
    B_level = B_level - P_flight_res(p_count)*ts; % Add energy gained during 'ts' time to battery level. Since P_flight_res is negative, subtract it. 
    t_rech_batt = t_rech_batt + ts; % Increment time. 
    if p_count == dusk_start_index
%         disp("Batteries did not recharge to last the next night cycle")
    end
    p_count = p_count + 1; % Increment p-counter
end

[~, ~, v] = Drag_Polar(aero, geom, S, W_t); % Find required speed to output. 

% Mass stuff
M_total = W_t/9.81;                 % Find total mass.
M_e     = (W_t-W_b-PT.W_s)/9.81;    % Find empty mass. 

% Finding times
t_night  = t_night/3600;                        % Convert time to hours
t_day    = t_day/3600;                          % Convert time to hours
t_batt = (PT.eta_batt*(E_b ...
    - Echem_dawn - Echem_dusk)/P_flight)/3600;  % Time of flight using only battery [hours].
t_rech_night = t_rech_night/3600;               % Convert time to hours
t_rech_batt = t_rech_batt/3600;                 % Convert time to hours

end


