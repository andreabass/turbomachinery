% %%NACA 65 SERIES BLADE DATA
% x_perc=[0 1.25 2.5 5 7.5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percentage x/c
% y_t_perc=[0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151 3.038 1.847 0.749 0.354 0.15];% half thickness t/c
% y_c_perc=[0 0.535 0.93 1.58 2.12 2.585 3.365 3.98 4.86 5.355 5.515 5.355 4.86 3.98 2.585 1.58 0];% camber line y/c
% dy_dx_perc=[Inf 0.3477 0.2915 0.2343 0.1999 0.1749 0.1381 0.1103 0.0675 0.0323 0 0.0323 0.0675 0.1103 0.1749 0.2343 -Inf];% derivative
% 

%% DOUBLE CIRCULAR ARC BLADE PROFILE
x_perc=[0 1.25 2.5 5 7.5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percentage x/c
y_t_perc=[0 0.217 0.334 0.534 0.778 0.985 1.369 1.711 2.275 2.678 2.9195 3 2.819 2.274 1.364 0.771 0.1];% half thickness t/c

epsi_h = 0.5 * tand(((beta_2_h_geo-beta_1_h_geo))/4);
R_perc_h = epsi_h/2 + 1/(8*epsi_h);
y_c_perc_h = epsi_h - R_perc_h + ( R_perc_h^2 - (x_perc./100 - 0.5).^2).^0.5;
dy_dx_perc_h = - (x_perc./100 - c_R_m/2).*(R_perc_h^2 - (x_perc./100 - 1/2).^2).^(-0.5);

epsi_m = 0.5 * tand(((beta_2_m_geo-beta_1_m_geo))/4);
R_perc_m = epsi_m/2 + 1/(8*epsi_m);
y_c_perc_m = epsi_m - R_perc_m + ( R_perc_m^2 - (x_perc./100 - 0.5).^2).^0.5;
dy_dx_perc_m = - (x_perc./100 - c_R_m/2).*(R_perc_m^2 - (x_perc./100 - 1/2).^2).^(-0.5);

epsi_t = 0.5 * tand(((beta_2_t_geo-beta_1_t_geo))/4);
R_perc_t = epsi_t/2 + 1/(8*epsi_t);
y_c_perc_t = epsi_t - R_perc_t + ( R_perc_t^2 - (x_perc./100 - 0.5).^2).^0.5;
dy_dx_perc_t = - (x_perc./100 - c_R_m/2).*(R_perc_t^2 - (x_perc./100 - 1/2).^2).^(-0.5);


% THICKNESS PROFILE
x_r = x_perc*c_R_m;
y_t_r = c_R_m*y_t_perc*0.8;


% MID-CHORD POINT
x_rAV = 50*c_R_m;
% cl_rm=4*pi*epsi_m/log(2);
% y_rAV = c_R_m*y_c_perc_m(15)*100*cl_rm;
y_rAV =0;
%% HUB 

% CAMBER LINE PROFILE
cl_rh= 4*pi*epsi_h/log(2);
y_rh = c_R_m*y_c_perc_h*100;
dy_dx_rh = dy_dx_perc_h*c_R_m;
eps_rh = atand(dy_dx_rh);

% SUCTION SIDE PROFILE (X-Y)
xss_rh = x_r-y_t_r.*sind(eps_rh);
yss_rh = y_rh+y_t_r.*cosd(eps_rh);

% PRESSURE SIDE PROFILE (X-Y)
xps_rh = x_r+y_t_r.*sind(eps_rh);
yps_rh = y_rh-y_t_r.*cosd(eps_rh);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_rh = -(xss_rh*sind(-gamma_roth)+yss_rh*cosd(-gamma_roth));
AX_coord_ss_rh = xss_rh*cosd(-gamma_roth)-yss_rh*sind(-gamma_roth);

T_coord_ps_rh = -(xps_rh*sind(-gamma_roth)+yps_rh*cosd(-gamma_roth));
AX_coord_ps_rh = xps_rh*cosd(-gamma_roth)-yps_rh*sind(-gamma_roth);

T_coord_AVrh = -(x_rAV*sind(-gamma_roth)+y_rAV*cosd(-gamma_roth));
AX_coord_AVrh = x_rAV*cosd(-gamma_roth)-y_rAV*sind(-gamma_roth);


%% MID 

% CAMBER LINE PROFILE
cl_rm=4*pi*epsi_m/log(2);
y_rm = c_R_m*y_c_perc_m*100;
dy_dx_rm = dy_dx_perc_m*cl_rm;
eps_rm = atand(dy_dx_rm);

% SUCTION SIDE PROFILE (X-Y)
xss_rm = x_r-y_t_r.*sind(eps_rm);
yss_rm = y_rm+y_t_r.*cosd(eps_rm);

% PRESSURE SIDE PROFILE (X-Y)
xps_rm = x_r+y_t_r.*sind(eps_rm);
yps_rm = y_rm-y_t_r.*cosd(eps_rm);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_rm = -(xss_rm*sind(-gamma_rotm)+yss_rm*cosd(-gamma_rotm));
AX_coord_ss_rm = xss_rm*cosd(-gamma_rotm)-yss_rm*sind(-gamma_rotm);

T_coord_ps_rm = -(xps_rm*sind(-gamma_rotm)+yps_rm*cosd(-gamma_rotm));
AX_coord_ps_rm = xps_rm*cosd(-gamma_rotm)-yps_rm*sind(-gamma_rotm); 

T_coord_AVrm = -(x_rAV*sind(-gamma_rotm)+y_rAV*cosd(-gamma_rotm));
AX_coord_AVrm = x_rAV*cosd(-gamma_rotm)-y_rAV*sind(-gamma_rotm);

%% TIP 

% CAMBER LINE PROFILE
cl_rt=4*pi*epsi_t/log(2);
y_rt = c_R_m*y_c_perc_t*100;
dy_dx_rt = dy_dx_perc_t*cl_rt;
eps_rt = atand(dy_dx_rt);

% SUCTION SIDE PROFILE (X-Y)
xss_rt = x_r-y_t_r.*sind(eps_rt);
yss_rt = y_rt+y_t_r.*cosd(eps_rt);

% PRESSURE SIDE PROFILE (X-Y)
xps_rt = x_r+y_t_r.*sind(eps_rt);
yps_rt = y_rt-y_t_r.*cosd(eps_rt);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_rt = -(xss_rt*sind(-gamma_rott)+yss_rt*cosd(-gamma_rott));
AX_coord_ss_rt = xss_rt*cosd(-gamma_rott)-yss_rt*sind(-gamma_rott);

T_coord_ps_rt = -(xps_rt*sind(-gamma_rott)+yps_rt*cosd(-gamma_rott));
AX_coord_ps_rt = xps_rt*cosd(-gamma_rott)-yps_rt*sind(-gamma_rott); 

T_coord_AVrt = -(x_rAV*sind(-gamma_rott)+y_rAV*cosd(-gamma_rott));
AX_coord_AVrt = x_rAV*cosd(-gamma_rott)-y_rAV*sind(-gamma_rott);

%% TRANSLATION OF HUB AND TIP PROFILES

% HUB TRANSLATION
DeltaT_AVrh = T_coord_AVrm - T_coord_AVrh;
DeltaAX_AVrh = AX_coord_AVrm - AX_coord_AVrh;

T_coord_ss_rh_trans = -(xss_rh*sind(-gamma_roth)+yss_rh*cosd(-gamma_roth))+DeltaT_AVrh;
AX_coord_ss_rh_trans = xss_rh*cosd(-gamma_roth)-yss_rh*sind(-gamma_roth)+DeltaAX_AVrh;

T_coord_ps_rh_trans = -(xps_rh*sind(-gamma_roth)+yps_rh*cosd(-gamma_roth))+DeltaT_AVrh;
AX_coord_ps_rh_trans = xps_rh*cosd(-gamma_roth)-yps_rh*sind(-gamma_roth)+DeltaAX_AVrh;

% TIP TRANSLATION
DeltaT_AVrt = T_coord_AVrm - T_coord_AVrt;
DeltaAX_AVrt = AX_coord_AVrm - AX_coord_AVrt;

T_coord_ss_rt_trans = -(xss_rt*sind(-gamma_rott)+yss_rt*cosd(-gamma_rott))+DeltaT_AVrt;
AX_coord_ss_rt_trans = xss_rt*cosd(-gamma_rott)-yss_rt*sind(-gamma_rott)+DeltaAX_AVrt;

T_coord_ps_rt_trans = -(xps_rt*sind(-gamma_rott)+yps_rt*cosd(-gamma_rott))+DeltaT_AVrt;
AX_coord_ps_rt_trans = xps_rt*cosd(-gamma_rott)-yps_rt*sind(-gamma_rott)+DeltaAX_AVrt;