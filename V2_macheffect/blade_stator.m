%% NACA 65 SERIES BLADE DATA
x_perc=[0 1.25 2.5 5 7.5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percentage x/c
y_t_perc=[0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151 3.038 1.847 0.749 0.354 0.15];% half thickness t/c
y_c_perc=[0 0.535 0.93 1.58 2.12 2.585 3.365 3.98 4.86 5.355 5.515 5.355 4.86 3.98 2.585 1.58 0];% camber line y/c
dy_dx_perc=[Inf 0.3477 0.2915 0.2343 0.1999 0.1749 0.1381 0.1103 0.0675 0.0323 0 0.0323 0.0675 0.1103 0.1749 0.2343 -Inf];% derivative

% THICKNESS PROFILE
x_s = x_perc*c_S_m;
y_t_s = c_S_m*y_t_perc*0.8; % Half-thickness

% MID-CHORD POINT
x_sAV = 50*c_S_m;
y_sAV = 0;

%% HUB 

% CAMBER LINE PROFILE
cl_sh=teta_hub_stat/25;         % Lift coefficient
y_sh = c_S_m*y_c_perc*cl_sh;      % Circular arc camber line
dy_dx_sh = dy_dx_perc*cl_sh;    % Circular arc camber line slope
eps_sh = atand(dy_dx_sh);       % Blade angles

% SUCTION SIDE PROFILE (X-Y)
xss_sh = x_s-y_t_s.*sind(eps_sh); 
yss_sh = y_sh+y_t_s.*cosd(eps_sh);

% PRESSURE SIDE PROFILE (X-Y)
xps_sh = x_s+y_t_s.*sind(eps_sh);
yps_sh = y_sh-y_t_s.*cosd(eps_sh);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_sh = -(xss_sh*sind(gamma_stath)+yss_sh*cosd(gamma_stath));
AX_coord_ss_sh = xss_sh*cosd(gamma_stath)-yss_sh*sind(gamma_stath);

T_coord_ps_sh = -(xps_sh*sind(gamma_stath)+yps_sh*cosd(gamma_stath));
AX_coord_ps_sh = xps_sh*cosd(gamma_stath)-yps_sh*sind(gamma_stath);

T_coord_AVsh = -(x_sAV*sind(gamma_stath)+y_sAV*cosd(gamma_stath));
AX_coord_AVsh = x_sAV*cosd(gamma_stath)-y_sAV*sind(gamma_stath);

%% MID 

% CAMBER LINE PROFILE
cl_sm=teta_mid_stat/25;
y_sm = c_S_m*y_c_perc*cl_sm;
dy_dx_sm = dy_dx_perc*cl_sm;
eps_sm = atand(dy_dx_sm);

% SUCTION SIDE PROFILE (X-Y)
xss_sm = x_s-y_t_s.*sind(eps_sm);
yss_sm = y_sm+y_t_s.*cosd(eps_sm);

% PRESSURE SIDE PROFILE (X-Y)
xps_sm = x_s+y_t_s.*sind(eps_sm);
yps_sm = y_sm-y_t_s.*cosd(eps_sm);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_sm = -(xss_sm*sind(gamma_statm)+yss_sm*cosd(gamma_statm));
AX_coord_ss_sm = xss_sm*cosd(gamma_statm)-yss_sm*sind(gamma_statm);

T_coord_ps_sm = -(xps_sm*sind(gamma_statm)+yps_sm*cosd(gamma_statm));
AX_coord_ps_sm = xps_sm*cosd(gamma_statm)-yps_sm*sind(gamma_statm); 

T_coord_AVsm = -(x_sAV*sind(gamma_statm)+y_sAV*cosd(gamma_statm));
AX_coord_AVsm = x_sAV*cosd(gamma_statm)-y_sAV*sind(gamma_statm);

%% TIP 

% CAMBER LINE PROFILE
cl_st=teta_tip_stat/25;
y_st = c_S_m*y_c_perc*cl_st;
dy_dx_st = dy_dx_perc*cl_st;
eps_st = atand(dy_dx_st);

% SUCTION SIDE PROFILE (X-Y)
xss_st = x_s-y_t_s.*sind(eps_st);
yss_st = y_st+y_t_s.*cosd(eps_st);

% PRESSURE SIDE PROFILE (X-Y)
xps_st = x_s+y_t_s.*sind(eps_st);
yps_st = y_st-y_t_s.*cosd(eps_st);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_st = -(xss_st*sind(gamma_statt)+yss_st*cosd(gamma_statt));
AX_coord_ss_st = xss_st*cosd(gamma_statt)-yss_st*sind(gamma_statt);

T_coord_ps_st = -(xps_st*sind(gamma_statt)+yps_st*cosd(gamma_statt));
AX_coord_ps_st = xps_st*cosd(gamma_statt)-yps_st*sind(gamma_statt); 

T_coord_AVst = -(x_sAV*sind(gamma_statt)+y_sAV*cosd(gamma_statt));
AX_coord_AVst = x_sAV*cosd(gamma_statt)-y_sAV*sind(gamma_statt);

%% TRANSLATION OF HUB AND TIP PROFILES

% HUB TRANSLATION
DeltaT_AVsh = T_coord_AVsm - T_coord_AVsh;
DeltaAX_AVsh = AX_coord_AVsm - AX_coord_AVsh;

T_coord_ss_sh_trans = T_coord_ss_sh+DeltaT_AVsh;
AX_coord_ss_sh_trans = AX_coord_ss_sh+DeltaAX_AVsh;

T_coord_ps_sh_trans = T_coord_ps_sh+DeltaT_AVsh;
AX_coord_ps_sh_trans = AX_coord_ps_sh+DeltaAX_AVsh;

% TIP TRANSLATION
DeltaT_AVst = T_coord_AVsm - T_coord_AVst;
DeltaAX_AVst = AX_coord_AVsm - AX_coord_AVst;

T_coord_ss_st_trans = T_coord_ss_st+DeltaT_AVst;
AX_coord_ss_st_trans = AX_coord_ss_st+DeltaAX_AVst;

T_coord_ps_st_trans = T_coord_ps_st+DeltaT_AVst;
AX_coord_ps_st_trans = AX_coord_ps_st+DeltaAX_AVst;

%% MIRRORING

T_coord_ss_st_trans = -T_coord_ss_st_trans;
T_coord_ps_st_trans = -T_coord_ps_st_trans;

T_coord_ss_sh_trans = -T_coord_ss_sh_trans;
T_coord_ps_sh_trans = -T_coord_ps_sh_trans;

T_coord_ss_sm = -T_coord_ss_sm;
T_coord_ps_sm = -T_coord_ps_sm;

