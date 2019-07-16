clc
clear
close all

%% OPTIONS
options

%% INPUT VALUES

% HUB/TIP RATIO
lambda = 0.62;

% ROTATIONAL SPEED
n = 5300; %rpm

% ROTOR INLET RELATIVE MACH NUMBER @ TIP
M_R1_t = 0.99;

% STATOR FLOW DEFLECTION
V2Tm_V1Tm = 0.5; % Put 0 for axial velocity at stator outlet

% SOLIDITIES 
sigma_R_m_design = 1.3;
sigma_S_m_design = 1.3;

% CHORDS
c_IGV = 0.04;
c_R_design = 0.06;
c_S_design = 0.05;

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

%% TRIANGLES

% MACHINE
tvd = figure;
subplot(2,2,1)
velt(V_0,0,0)
title('SECTION 0')
subplot(2,2,2)
velt(V_1_h,W_1_h,U_h,'k--')
velt(V_1_m,W_1_m,U_m,'k-.')
velt(V_1_t,W_1_t,U_t)
title('SECTION 1')
subplot(2,2,3)
velt(V_2_h,W_2_h,U_h,'k--')
velt(V_2_m,W_2_m,U_m,'k-.')
velt(V_2_t,W_2_t,U_t)
title('SECTION 2')
subplot(2,2,4)
velt(V_3A(end),V_3T_h,0,'k--')
velt(V_3A(end),V_3T_m,0,'k-.')
velt(V_3A(end),V_3T_t,0)
title('SECTION 3')

% ROTOR
tvdrot = figure;
subplot(3,1,3)
velt(V_1_h,W_1_h,U_h)
velt(V_2_h,W_2_h,U_h,'k--')
title('ROTOR HUB')
subplot(3,1,2)
velt(V_1_m,W_1_m,U_m)
velt(V_2_m,W_2_m,U_m,'k--')
title('ROTOR MID')
subplot(3,1,1)
velt(V_1_t,W_1_t,U_t)
velt(V_2_t,W_2_t,U_t,'k--')
title('ROTOR TIP')





