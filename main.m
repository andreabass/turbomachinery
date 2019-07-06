clc
clear
close all

%% OPTIONS
options

%% INPUT VALUES
D_t = 1; %m
lambda = 0.569;
n = 9300; %rpm
omega = n * 2 * pi / 60; %rad/s

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

geo_rotor
geo_stator
blade_rotor
blade_stator



