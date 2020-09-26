function [e_tot, p_vec, t_day, t_night, day_vec, ts] = Solar_Array_Model(PT)
%Models the solar energy captured by the solar array. e_tot is J/m^2.
%P_vec is a vector of the power output over a 24 hr cycle [W/m^2].
%t_day and t_night are scalars of the total duration of the day/night [s].
%day_vec is a vector of times throughout the day [seconds]. 
%ts is the TimeStep between data points in the textfile in seconds. Data based off
%NASA CR 3699 research paper.
%Luke Bughman, 9/13/18

%% read in sun data
FID = fopen(strcat(PT.SunData,'.txt'),'r');          %open the file

%SunData = textscan(FID,'%{hh:mm}T%f%f');             %read in file using textscan - UNCOMMENT THIS LINE ONLY IF RUNNING R2018a OR LATER
%clock_time = SunData{:,1};                           %UNCOMMENT THIS LINE ONLY IF RUNNING R2018a OR LATER

SunData = textscan(FID,'%f%c%f%f%f');                   %uncomment if using pre-R2018a versions
clock_time = hours(SunData{:,1})+minutes(SunData{:,3});   %uncomment if using pre-R2018a versions
theta_temp = SunData{1,4};                           %theta is sun elevation angle

%time and elevation angle after sunrise
time = clock_time(theta_temp >= 0);                  %time after sunrise
theta = theta_temp(theta_temp >= 0);                 %elevation after sunrise

%% add data for atmospeheric attenuation at sea level
theta_c_a = [0 5 10 15 20 30 45 90];                    %This data is from NASA CR 3699
c_a_interp = [0 .1 .2 .31 .38 .48 .56 .65];

I_nom = 226.56;            %[W/m^2] nominal solar power, no atmospheric attenuation. July 5th. 
                           %Old Value: 1309. New value already takes into
                           %account the %20 efficiency of the panels, which
                           %is why it is so much lower. 

%% Calculate power recieved by solar array throught the day
for i = 1:length(time)
    phi = abs(theta(i)-PT.array_norm);      %Phi is the angle between array normal and sun [degrees]
    c_a = interp1(theta_c_a,c_a_interp,theta(i));     %Interpolate to find atmospheric attenuation
    I_act(i) = c_a*I_nom*cosd(phi);            %actual power recived by solar array [w/m^2]
    
    if phi >= 90
        I_act(i) = NaN;                %sun is past solar cell - cannot receive power
    end
end

%% use solar array efficiency to calculate power out of solar array;
p_vec = PT.eta_mppt*PT.eta_solar.*I_act;                %Output of solar array AFTER MPPT [W/m^2]

%% Integrate over the course of the day to find total energy per meter squared
%find timestep of data
ts = seconds(time(2))-seconds(time(1));           %[s]
day_vec = seconds(time)-seconds(time(1));                       %elapsed seconds since sunrise
e_tot = PT.eta_mppt*max(cumtrapz(day_vec,p_vec));                            %integrate to find total energy per meter squared throughtout the day [J/m^2]
e_tot = e_tot*PT.A;         % Find total energy per solar panel [J/panel].
p_vec = p_vec*PT.A;         % Find power vector per solar panel [W/panel].

% for i = 1:length(day_vec)
%     delta_e(i) = P_vec(i)*ts;   %change in energy with each timestep [J]
%     e_tot_f(i) = sum(delta_e(1:i));
% end
% 
% e_tot = e_tot_f(end);

%find time during day and night for other functions
t_day = max(day_vec);
t_night = 24*3600-t_day;

%% close all files
fclose('all');
end

