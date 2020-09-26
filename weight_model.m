function [W_e] = weight_model(~,b,~,~,geom,aero,~)
%Preliminary Weight Estimates for Project Mobius
%all Weights in Newtons
%work in progress - will be updated as components are refined and test
%sections built

%% Constant Weight Parameters

W_const = 4.16925;  % Taken from prototype 1 parts sheet.[N]

% Mobius's stuff.
% W_servos = .94;                 %4x 24 gram servos
% W_MPPT = 3;                     %Pat's new quad MPPT
% W_cam = .4;                    %camera weight
% W_avionics = 1.3;                %from drive
% W_spinner = .41;                %currently used on mark 3

% %Sum constant weights
% W_const = W_servos+W_MPPT+W_cam+W_avionics+W_spinner;


%% non constant weight parameters
W_wire = 9.81*.06*b;    %wire weight is a function of wingspan 
W_prop = 0.2;           % Pure guess (~20 g)
% W_prop = .027*W_to;     %based on Mk3 Prototype prop
% W_motor = .0046*W_to*design.PWto;   %based on Mk3 motor
W_motor = 0;            % Already accounted for in constant weights. 
% W_skid = .007*W_to;   % will we even have a skid??
W_skid = 0;             % Assuming no skid. 

%sum to find variable weights
W_var = W_wire+W_prop+W_motor+W_skid;

%% Structural Weight Estimation
MGC = b/aero.AR;                    %mean geometric chord
b_joint = .1*b;                     %for joined wing truss configuration - assuming rear wing is half the span of main wing
t = geom.t_c*MGC;                   %wing thickness

%balsa properties
rho_balsa = 160*9.81;        %[N/m^3]
t_balsa = 0.003175;                  %[m] [1/8" balsa];

%size ribs
n_ribs = round((b+b_joint)/(MGC*.5));    %calculate number of balsa ribs - spaced .5 chordlengths apart
A_rib = .6*t*MGC;                   %approximation for rib area
V_ribs = n_ribs*A_rib*t_balsa;       %[m^3]
W_ribs = 1.2*rho_balsa*V_ribs;     %[N] add 20% for epoxy weight           

%size the shear web
A_SW = (b_joint+b)*t;              %[m^2]
V_SW = A_SW*t_balsa;               %[m^3]
W_SW = 1.2*rho_balsa*V_SW;         %[N] add 20% for epoxy weight

%size fiberglass wing skin
L_skin = 1.2*MGC;                  %[m] section length of skin - won't wrap fully around wing
A_skin = L_skin*(b_joint+b);       %[m^2]
W_skin = 4.11*A_skin;               %[N] - based off test layups

%size monokote 
L_mono = .8*MGC;                    %[m] section length of monokote (on bottom surface of airfoil)
A_mono = L_mono*(b_joint+b);
W_mono = .6*A_mono;                %[N]

%size fuselage - based off guesses
L_fuse = .5*b;                       %[m]
W_fuse = 3*L_fuse;                   %[N]

%Size Vertical - assuming weight is proportional to wingspan
W_vert = .5*b;                       %[N]


% % size winglets for flying wing
% if aero.JW == 0
%     W_winglets = .77*S;
%     
%     %spar weight derivation is separate - constant in front of equation
%     %will be adjusted after a prototype wing section is created
%     W_spar = (6.72*10^-5)*(loads.n*W_to*b*aero.AR)/geom.t_c;     
% else
%     W_winglets = 0;
% end
% 
% W_struct = W_le+W_te+W_ribs+W_winglets+W_spar;

%W_struct = 12.82*S;           %from composite wing section

%W_struct = 10*S;           %from composite wing section

W_struct = W_ribs+W_SW+W_skin+W_mono+W_fuse+W_vert;     %[N]

%% Calculate Empty Weight and Empty Weight fraction
W_e = W_const+W_var+W_struct;
% ew_frac = W_e/W_to;

end

