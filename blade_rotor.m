% ROTOR
mod_gamma_roth = abs(beta_1_h)-abs(teta_hub)/2-(i_opt_roth);
gamma_roth =mod_gamma_roth;

mod_gamma_rotm = abs(beta_1_m)-abs(teta_mid)/2-(i_opt_rotm);
gamma_rotm =mod_gamma_rotm;

mod_gamma_rott = abs(beta_1_t)-abs(teta_tip)/2-(i_opt_rott);
gamma_rott =mod_gamma_rott;

%-------------------------------------------------------------------------
%                           NACA 65

x_perc=[0 1.25 2.5 5 7.5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percenage x/c
y_t_perc=[0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151 3.038 1.847 0.749 0.354 0.15];% half thickness t/c
y_c_perc=[0 0.535 0.93 1.58 2.12 2.585 3.365 3.98 4.86 5.355 5.515 5.355 4.86 3.98 2.585 1.58 0];% camber line y/c
dy_dx_perc=[Inf 0.3477 0.2915 0.2343 0.1999 0.1749 0.1381 0.1103 0.0675 0.0323 0 0.0323 0.0675 0.1103 0.1749 0.2343 -Inf];% derivative

%------------------------------------------------------------------------

c_r=c_R_m;

x_r = x_perc*c_r;
y_t_r = c_r*y_t_perc*0.8;
x_rAV = 50*c_r;
y_rAV = 0;

% HUB 

cl_rh=teta_hub/25;

% equations

y_rh = c_r*y_c_perc*cl_rh;
dy_dx_rh = dy_dx_perc*cl_rh;

eps_rh = atand(dy_dx_rh);

xss_rh = x_r-y_t_r.*sind(eps_rh);
yss_rh = y_rh+y_t_r.*cosd(eps_rh);

xps_rh = x_r+y_t_r.*sind(eps_rh);
yps_rh = y_rh-y_t_r.*cosd(eps_rh);

% Passage from x-y coordinates to T-ax

T_coord_ss_rh = -(xss_rh*sind(gamma_roth)+yss_rh*cosd(gamma_roth));
AX_coord_ss_rh = xss_rh*cosd(gamma_roth)-yss_rh*sind(gamma_roth);

T_coord_ps_rh = -(xps_rh*sind(gamma_roth)+yps_rh*cosd(gamma_roth));
AX_coord_ps_rh = xps_rh*cosd(gamma_roth)-yps_rh*sind(gamma_roth);

T_coord_AVrh = -(x_rAV*sind(gamma_roth)+y_rAV*cosd(gamma_roth));
AX_coord_AVrh = x_rAV*cosd(gamma_roth)-y_rAV*sind(gamma_roth);


% MID 

cl_rm=teta_mid/25;

% equations

y_rm = c_r*y_c_perc*cl_rm;
dy_dx_rm = dy_dx_perc*cl_rm;

eps_rm = atand(dy_dx_rm);

xss_rm = x_r-y_t_r.*sind(eps_rm);
yss_rm = y_rm+y_t_r.*cosd(eps_rm);

xps_rm = x_r+y_t_r.*sind(eps_rm);
yps_rm = y_rm-y_t_r.*cosd(eps_rm);

% Passage from x-y coordinates to T-ax

T_coord_ss_rm = -(xss_rm*sind(gamma_rotm)+yss_rm*cosd(gamma_rotm));
AX_coord_ss_rm = xss_rm*cosd(gamma_rotm)-yss_rm*sind(gamma_rotm);

T_coord_ps_rm = -(xps_rm*sind(gamma_rotm)+yps_rm*cosd(gamma_rotm));
AX_coord_ps_rm = xps_rm*cosd(gamma_rotm)-yps_rm*sind(gamma_rotm); 

T_coord_AVrm = -(x_rAV*sind(gamma_rotm)+y_rAV*cosd(gamma_rotm));
AX_coord_AVrm = x_rAV*cosd(gamma_rotm)-y_rAV*sind(gamma_rotm);


% TIP 

cl_rt=teta_tip/25;

% equations

y_rt = c_r*y_c_perc*cl_rt;
dy_dx_rt = dy_dx_perc*cl_rt;

eps_rt = atand(dy_dx_rt);

xss_rt = x_r-y_t_r.*sind(eps_rt);
yss_rt = y_rt+y_t_r.*cosd(eps_rt);

xps_rt = x_r+y_t_r.*sind(eps_rt);
yps_rt = y_rt-y_t_r.*cosd(eps_rt);

% Passage from x-y coordinates to T-ax

T_coord_ss_rt = -(xss_rt*sind(gamma_rott)+yss_rt*cosd(gamma_rott));
AX_coord_ss_rt = xss_rt*cosd(gamma_rott)-yss_rt*sind(gamma_rott);

T_coord_ps_rt = -(xps_rt*sind(gamma_rott)+yps_rt*cosd(gamma_rott));
AX_coord_ps_rt = xps_rt*cosd(gamma_rott)-yps_rt*sind(gamma_rott); 

T_coord_AVrt = -(x_rAV*sind(gamma_rott)+y_rAV*cosd(gamma_rott));
AX_coord_AVrt = x_rAV*cosd(gamma_rott)-y_rAV*sind(gamma_rott);

% HUB
DeltaT_AVrh = T_coord_AVrm - T_coord_AVrh;
DeltaAX_AVrh = AX_coord_AVrm - AX_coord_AVrh;

T_coord_ss_rh_trans = -(xss_rh*sind(gamma_roth)+yss_rh*cosd(gamma_roth))+DeltaT_AVrh;
AX_coord_ss_rh_trans = xss_rh*cosd(gamma_roth)-yss_rh*sind(gamma_roth)+DeltaAX_AVrh;

T_coord_ps_rh_trans = -(xps_rh*sind(gamma_roth)+yps_rh*cosd(gamma_roth))+DeltaT_AVrh;
AX_coord_ps_rh_trans = xps_rh*cosd(gamma_roth)-yps_rh*sind(gamma_roth)+DeltaAX_AVrh;

% TIP
DeltaT_AVrt = T_coord_AVrm - T_coord_AVrt;
DeltaAX_AVrt = AX_coord_AVrm - AX_coord_AVrt;

T_coord_ss_rt_trans = -(xss_rt*sind(gamma_rott)+yss_rt*cosd(gamma_rott))+DeltaT_AVrt;
AX_coord_ss_rt_trans = xss_rt*cosd(gamma_rott)-yss_rt*sind(gamma_rott)+DeltaAX_AVrt;

T_coord_ps_rt_trans = -(xps_rt*sind(gamma_rott)+yps_rt*cosd(gamma_rott))+DeltaT_AVrt;
AX_coord_ps_rt_trans = xps_rt*cosd(gamma_rott)-yps_rt*sind(gamma_rott)+DeltaAX_AVrt;

figure
% hubr = plot(T_coord_ss_rh_trans,-AX_coord_ss_rh_trans,'green',T_coord_ps_rh_trans,-AX_coord_ps_rh_trans,'green')
hubrss = plot(T_coord_ss_rh_trans,-AX_coord_ss_rh_trans,'green');
hold on
hubrps = plot(T_coord_ps_rh_trans,-AX_coord_ps_rh_trans,'green');
midrss = plot(T_coord_ss_rm,-AX_coord_ss_rm,'red');
midrps = plot(T_coord_ps_rm,-AX_coord_ps_rm,'red');
tiprss = plot(T_coord_ss_rt_trans,-AX_coord_ss_rt_trans,'black');
tiprps = plot(T_coord_ps_rt_trans,-AX_coord_ps_rt_trans,'black');
title('ROTOR BLADE (NACA 65)')
legend([hubrss midrss tiprss],{'HUB','MID','TIP'})
hold off

axis([-3.25 0.5 -2.95 0.8])

pbaspect([1 1 1])