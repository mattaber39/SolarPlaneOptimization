function [results] = Constant_Weight_Sizing(design, aero, PT, sys, geom, loads,vis)
%constant weight sizing for solar power aircraft
%this Iterative solver first finds the energy required to fly through the
%night. The batteries are then sized based off this. Then the energy required to
%fly throught the day is added to e_night, and the solar array is sized
%based off that requirement

%% extract variables
WtoS = design.WtoS;                %design wing loading [N/m^2]

%% Load in solar array model
PT.array_norm = 90;
[e_tot, p_vec, PT.t_day, PT.t_night,~,ts] = Solar_Array_Model(PT);

PT.array_norm = 0;
[e_tot_vert, p_vec_vert] = Solar_Array_Model(PT);

%e_tot is total energy per square meter from solar array [J/m^2]
%p_vec is a vector of array power throughout the day (for plotting purpose)
%[W/m^2]
%t_day and t_night are length of day and night in seconds

%% iteravely solve to find Takeoff Weight
Wto_old = 20;                %initial guess for takeoff weight [N] 
S_solar_old = 1;
done = false;                %logical used for while loop
tol  = .01;                  %tolerance on solver [N]
it = 0;                      %start iteration counter
max_it = 500;                %maximum number of iterations

while ~done 
    it = it+1;
    
    %find old wing area and span
    S = Wto_old/WtoS;
    b = sqrt(aero.AR*S);
    
    %size battery to fly throught the night (all efficiencies accounted for
    %in respective functions)
    [W_batt, E_batt, V_batt] = size_batt(Wto_old,S,aero,geom,PT,sys,p_vec,S_solar_old,ts); 
    
    %find energy required during the day (excluding battery charging)
    [P_thrust] = thrust_power_required(aero,geom,Wto_old,S);
    E_day = (P_thrust+sys.P)*PT.t_day;
    E_night = (P_thrust+sys.P)*PT.t_night/(PT.eta_cc*PT.batt_range*PT.eta_batt); 
    
    %find area and weight of solar cells
    E_tot = E_day+1.07*E_night;                   %total energy required FROM THE SOLAR CELLS [J];
    S_solar = E_tot/((1-PT.array_vert)*e_tot+PT.array_vert*e_tot_vert/2);                      %required area of solar collector [m^2]
%     N_solar = 
    W_solar = S_solar*PT.sigma;                 %[N]
    
    %use weight model to calculate empty weight
    [W_e] = weight_model(Wto_old,b,design,S,geom,aero,loads);
    
    %calculate new takeoff weight
    Wto = sys.W_p+W_solar+W_batt+W_e;
    
    %check if done
    done = abs(Wto_old-Wto) <= tol;
    
    %update guess
    Wto_old = Wto;     
    S_solar_old = S_solar;
    
    if it >= max_it
        break
    end
end

%% Store results in results structure
[results.C_L_cruise,C_D_cruise, results.V_min_pwr] = Drag_Polar(aero,geom,S,Wto);
results.LoD_cruise = results.C_L_cruise/C_D_cruise;
results.W_solar = W_solar;
results.S = S;
results.S_solar = S_solar;
results.Wto = Wto;
results.Wes = W_batt;
results.V_batt = V_batt;
results.solar_frac = results.S_solar/results.S;
% results.EW_frac = EW_frac;
results.E_batt = E_batt/3600;                    %total energy stored batteries [W-h]
results.b = b;

%% Plot Battery capacity
if vis == 1

figure()
hold on
plot((0:ts:PT.t_day)./60,(1-PT.array_vert)*S_solar*p_vec)
plot((0:ts:PT.t_day)./60,PT.array_vert*S_solar*p_vec_vert)
plot((0:ts:PT.t_day)./60,(1-PT.array_vert)*S_solar*p_vec+PT.array_vert*S_solar*p_vec_vert)
xlabel('Minutes after Sunrise')
ylabel('Power from Array [W]')
grid on

P_thrust = thrust_power_required(aero,geom,Wto,S);
    
t_vec = 0:ts:24*3600;                          %time vector [s]
P_vec = (1-PT.array_vert)*S_solar*p_vec+PT.array_vert*S_solar*p_vec_vert;                        %power produced by solar cells (MPPT efficiency already included)
P_const = -(P_thrust+sys.P);

for i = 1:length(t_vec)
    E_cons(i) = t_vec(i).*P_const;
    
    if i <= length(P_vec)
        E_solar(i) = ts*trapz(P_vec(1:i));
        E_batt_cycle(i) = (E_solar(i)+E_cons(i))*PT.eta_cc;
        chrg_end = E_batt_cycle(i);
    else
        E_solar(i) = NaN;
        E_batt_cycle(i) = chrg_end+(t_vec(i)-PT.t_day).*P_const;
    end
end
const = -min(E_batt_cycle);

figure()
hold on
plot(t_vec./3600,E_solar./3600)
plot(t_vec./3600, E_cons./3600)
plot(t_vec./3600, (E_batt_cycle+const)./3600)
xlabel('Hours after Sunrise')
ylabel('Energy [W-h]')
legend('Solar Power (Energy in)','Power to Motor (Energy Out)','Battery Storage')
grid on

end
