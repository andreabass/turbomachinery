clc
clear
close all

%% OPTIONS
options

%% INPUT VALUES

% 1
D_t = 1; %m
% 2
lambda = 0.55;
% 3
n = 10000; %rpm
% 4
c_IGV = 0.04;
% 5
t_over_s_h = 0.02; % Relative trailing edge thickness
% 6
th_c = 0.08; % Percentage thickness (same for IGV, rotor and stator blades
% 7
% V_3T_m
% 8 
sigma_R_m_design = 1;
% 9
sigma_S_m_design = 1.5;
% 10
c_R_design = 0.05;
% 11
c_S_design = 0.05;

%% PROBLEM DATA
data

%% DESIGN PROBLEM (PERFORMANCES)
eta_TT = 0.8;
eta_TT = [2*eta_TT eta_TT];

while abs((eta_TT(end)-eta_TT(end-1))/eta_TT(end-1))>tol
    eta_TT(end-1) = eta_TT(end);

%%% [ p1: IGV ] %%%
% 0 DOF

p1_pre
p1_igvinlet
p1_rotin_m_t_h
    
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





