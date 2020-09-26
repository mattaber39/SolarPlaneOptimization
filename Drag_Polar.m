function [C_L_min_pwr, C_D_min_pwr, V_min_pwr] = Drag_Polar(aero,geom,S,W)
%This drag polar only applies for Steady Level Flight, and is a function of 
%the aircraft velocity. Takes into account
%Parasitic, RE-dependent viscous, Lift dependent Viscous, and Induced drag.
%must specify actual wing loading and area (not WtoS), therefore must be
%used with an iterative solver. Find C_L C_D and Velocity for minimum power
%required. Does not include trim drag yet.

%% setup calcs
%find atmostpheric condictions
[~, ~, ~, rho] = atmosisa(aero.alt);
mu = 1.9*10^-5;          %dynamic viscosity @ 40 C

%% create vector of velocities (to find min pwr required speed)
C_L_low = .2;                         %bottom of C_L range for vector
V_max = sqrt((2*W)/(rho*C_L_low*S));      %maximum speed aircraft could fly at
V_min = sqrt((2*W)/(rho*aero.C_L_max*S));  %minimum speed aircraft could fly at
V_vect = linspace(V_min,V_max);            %create vector of potential speeds aircraft could fly at

%% Iterate through speeds to find the minimum power required speed
done = false;
it = 1;
pwr = 0;

while ~done
    
   %update counter
   it = it+1;
   
   %if minimum power required speed hasn't been found, set results to NaN
   if it > length(V_vect)
       C_L_min_pwr = NaN;
       C_D_min_pwr = NaN;
       V_min_pwr = NaN;
       break
   end
   
    %% calculate Re number
    chord = sqrt(S/aero.AR)/2;         %SMC
    Re = (rho*V_vect(it)*chord)/mu;

    %find C_L that A/C is flying at in SLF
    C_L(it) = W/(.5*rho*V_vect(it)^2*S);

    %% calulcate drag terms
    C_DP = aero.C_DP;             %miscellaneous Parasitic drag term
    
    %Wing parasitic drag
    C_D_min = 1.0603/Re^.3645;    %for E174 airfoil
    C_L_min = .4-8*C_D_min;
    C_DP_wing = aero.k_dp*(C_L(it)-C_L_min)^2+C_D_min;      %parasitic drag of wing
    
    %Skin Friction Drag
    S_wet = 2.1*(S+2*aero.h*S);                               %wetted surface area of joined wing
    C_f = .455/(log(Re)^2.58);                                %skin coefficient of friction
    FF = 1+2*(geom.t_c)+60*(geom.t_c)^4;                      %Form factor for wing
    C_D_f = (C_f*FF*S_wet)/S;                                 %total skin friction drag

    %induced drag term
    C_D_ind = (C_L(it))^2/(pi*aero.e*aero.AR);

    %sum all components for total drag;
    C_D(it) = C_DP+C_DP_wing+C_D_f+C_D_ind;  
    
    %find pwr
    pwr(it) = C_L(it)^1.5/C_D(it);
    
    %check if done
    done = pwr(it) < pwr(it-1);
end

C_L_min_pwr = C_L(it-1);
C_D_min_pwr = C_D(it-1);
V_min_pwr = V_vect(it-1);

end

