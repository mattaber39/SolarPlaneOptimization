%% Function To Optimize
function [t_batt, t_extra, t_rech_night, t_rech_batt, t_day, t_night ...
    n_s, M_total, M_e, v, t_rech_night_extra] = function2Optimize(n_b, b, AR, aero, PT, geom)
% Still messy. Will clean up soon. 

aero.AR = AR;
S = b^2/aero.AR;            % Area [m2].
E_b = n_b*PT.e_b;           % Total energy of battery pack [J].
W_b = n_b*PT.w_b;           % Find weight of batteries [N].
W_e = weight_model...       % Find empty weight using Lukes' weight model [N].
    (0,b,0,0,geom,aero,0);

% Find the total number of solar cells able to fit on the wing. 
n_s = get_panels(b, aero.AR);  
PT.W_s = n_s*PT.wsps;

PT.W_s = n_s*PT.wsps;
W_t = W_e+W_b+PT.W_s;       % Find total estimated weight [N].
                            % PT.W_s = estimated weight of solar panels.

theta_climb = 5*pi/180;    
theta_fall  = -5*pi/180; 
theta_level = 0;

% Find thrust required. 
P_flight_level = thrust_power_required(aero,geom,W_t,S, theta_level);

% Find e_s, the energy each solar panel creates over the course of the day,
% in [W]. Find t_day, t_night (the length of day and night), in seconds. 
% Find the vector of power given by each solar panel throughout the day.
% Also find ts, the time step between each element of power. 
[e_s, p_vec, t_day, t_night, ~, ts] = Solar_Array_Model(PT);
p_s = e_s/t_day;

% Find power vector given by the whole solar array throughout the day. 
P_vec = n_s*p_vec;

% Find the average power throughout the day. 
P_s   = n_s*p_s;

% Find residual of power required to fly after subtracting power input by
% the solar array. 
P_flight_res = P_flight_level - P_vec;

Emech_night = P_flight_level*t_night;           % Thrust energy need to fly during night.
Echem_night = Emech_night/PT.eta_batt;    % Battery energy needed to fly during night.

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
Echem_dawn  = Emech_dawn/PT.eta_batt;     % Battery energy needed to fly during dawn.
Echem_dusk  = Emech_dusk/PT.eta_batt;     % Battery energy needed to fly during dusk.

if (Echem_dusk + Echem_night + Echem_dawn) < e_s
% Finding how long it takes to recharge for night. Checking if solar panels gather enough energy.
B_level = 0;            % Energy level of battery pack [J]
t_rech_night = 0;       % time to recharge for night only. 
p_count = t_dawn/ts+1;  % Start when solar array is giving excess power. 
while B_level < Echem_dusk + Echem_night + Echem_dawn && p_count <= dusk_start_index % While battery isn't charged for night and dusk hasn't arrived:
    B_level = B_level - P_flight_res(p_count)*ts;  % Add energy gained during 'ts' time to battery level. Since P_flight_res is negative, subtract it. 
    t_rech_night = t_rech_night + ts;  % Increment time. 
    p_count = p_count + 1;      % Increment p-counter
end

% Finding out how long it takes to fully recharge.
t_rech_batt  = t_rech_night;    % Start where we left off. 
while B_level < E_b && p_count <= dusk_start_index % While battery isn't full and dusk hasn't arrived:
    B_level = B_level - P_flight_res(p_count)*ts; % Add energy gained during 'ts' time to battery level. Since P_flight_res is negative, subtract it. 
    t_rech_batt = t_rech_batt + ts; % Increment time. 
    p_count = p_count + 1; % Increment p-counter
end

else
    t_rech_night = (Echem_dusk + Echem_night + Echem_dawn)/P_s;
    t_rech_batt  = E_b/P_s;
end

[~, ~, v] = Drag_Polar(aero, geom, S, W_t); % Find required speed to output. 

% Mass stuff
M_total = W_t/9.81;                 % Find total mass.
M_e     = (W_t-W_b-PT.W_s)/9.81;    % Find empty mass. 

% Finding times
t_night  = t_night/3600;                        % Convert time to hours
t_day    = t_day/3600;                          % Convert time to hours
t_batt = (PT.eta_batt*(E_b)/P_flight_level)/3600;  % Time of flight using only battery [hours].
t_rech_night = t_rech_night/3600;               % Convert time to hours
t_rech_batt = t_rech_batt/3600;                 % Convert time to hours
t_extra = t_batt - t_night - (t_dusk + t_dawn)/3600;
t_rech_night_extra = (t_day - (t_dusk + t_dawn)/3600) - t_rech_night;

end