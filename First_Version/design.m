clc
clear
close all

%% OPTIONS
options

%% INPUT VALUES

% 1
D_t = 1; %m
% 2
lambda = 0.50712;
% 3
n = 11000; %rpm
% 4
c_IGV = 0.04;
% 5
t_over_s_h = 0.02; % Relative trailing edge thickness
% 6
th_c = 0.08; % Percentage thickness (same for IGV, rotor and stator blades
% 7
V3Tm_V1Tm = 1; % Put 0 for axial velocity at stator outlet
% 8 
sigma_R_m_design = 1;
% 9
sigma_S_m_design = 1.7;
% 10
c_R_design = 0.05;
% 11
c_S_design = 0.05;

%% PROBLEM DATA
data

%% DESIGN PROBLEM (PERFORMANCES)
eta_TT = 0.87049;
eta_TT = [2*eta_TT eta_TT];

while abs((eta_TT(end)-eta_TT(end-1))/eta_TT(end-1))>tol
    eta_TT(end-1) = eta_TT(end);

%%% [ p1: IGV ] %%%
% 0 DOF

if DR==1
Xi_m=0.5;

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
geo_rotor
geo_stator
blade_IGV
blade_rotor
blade_stator
bladeplot
frontview


close all




