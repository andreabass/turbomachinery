clc
clear
close all

%% OPTIONS
options

%% INPUT VALUES

% HUB/TIP RATIO
lambda = 0.62;

% ROTATIONAL SPEED
n = 5250; %rpm

% ROTOR INLET RELATIVE MACH NUMBER @ TIP
M_R1_t = 0.99;

% STATOR FLOW DEFLECTION
V2Tm_V1Tm = 0.5; % Put 0 for axial velocity at stator outlet

% SOLIDITIES 
sigma_R_m_design = 1.3;
sigma_S_m_design = 1.3;

% CHORDS
c_IGV = 0.07;
c_R_design = 0.1;
c_S_design = 0.1;

% EXTERNAL DIAMETER
D_t = 1; %m

% OTHER GEOMETRICAL QUANTITIES
t_over_s_h = 0.02;  % Relative trailing edge thickness
th_c = 0.08;        % Percentage thickness (same for IGV, rotor and stator blades)

%% PROBLEM DATA
data

%% DESIGN PROBLEM (PERFORMANCES)
eta_TT = 0.95;
eta_TT = [2*eta_TT eta_TT];

while abs((eta_TT(end)-eta_TT(end-1))/eta_TT(end-1))>tol
    eta_TT(end-1) = eta_TT(end);

%%% [ p1: IGV ] %%%
% 0 DOF

if DR==1
Xi_m=0.6;

lambda = [2*lambda lambda];
while abs((lambda(end)-lambda(end-1))/lambda(end-1))>tol
    lambda(end-1) = lambda(end);
p1_pre_dr
p1_igvinlet_dr
p1_rotin_m_t_h_dr

end

        b = b(end);
        rho_0   = rho_0(end);
        rho_1_m = rho_1_m(end);
        rho_1_t = rho_1_t(end);
        rho_0   = rho_0(end);

else
    
p1_pre
p1_igvinlet
p1_rotin_m_t_h

end
        
%%% [ p2: ROTOR  ] %%%
% 1 DOF
% INPUT 1: eta_TT_m

p2_rotout_m_t_h

%%% [ p3: STATOR ] %%%
% 0 DOF

p3_statout_m_t_h

%%% [ UPDATE eta_TT,m ] %%%

updateval

end

pp_p1
pp_p2
pp_p3
check
results

%% DESIGN PROBLEM (GEOMETRY)

geo_IGV2
geo_rotor_final
geo_stator_final
blade_IGV
blade_rotor
blade_stator
bladeplot
frontview








