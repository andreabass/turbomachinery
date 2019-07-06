% STATOR
mod_gamma_stath = abs(alpha_2_h)-abs(teta_hub_stat)/2+abs(i_opt_stath);
gamma_stath =mod_gamma_stath;

mod_gamma_statm = abs(alpha_2_m)-abs(teta_mid_stat)/2+abs(i_opt_statm);
gamma_statm =mod_gamma_statm;

mod_gamma_statt = abs(alpha_2_t)-abs(teta_tip_stat)/2+abs(i_opt_statt);
gamma_statt =mod_gamma_statt;

c_s=c_S_m;

x_s = x_perc*c_s;
y_t_s = c_s*y_t_perc*0.8;
x_sAV = 50*c_s;
y_sAV = 0;

% HUB 

cl_sh=teta_hub_stat/25;

% equations

y_sh = c_s*y_c_perc*cl_sh;
dy_dx_sh = dy_dx_perc*cl_sh;

eps_sh = atand(dy_dx_sh);

xss_sh = x_s-y_t_s.*sind(eps_sh);
yss_sh = y_sh+y_t_s.*cosd(eps_sh);

xps_sh = x_s+y_t_s.*sind(eps_sh);
yps_sh = y_sh-y_t_s.*cosd(eps_sh);

% Passage from x-y coordinates to T-ax

T_coord_ss_sh = -(xss_sh*sind(gamma_stath)+yss_sh*cosd(gamma_stath));
AX_coord_ss_sh = xss_sh*cosd(gamma_stath)-yss_sh*sind(gamma_stath);

T_coord_ps_sh = -(xps_sh*sind(gamma_stath)+yps_sh*cosd(gamma_stath));
AX_coord_ps_sh = xps_sh*cosd(gamma_stath)-yps_sh*sind(gamma_stath);

T_coord_AVsh = -(x_sAV*sind(gamma_stath)+y_sAV*cosd(gamma_stath));
AX_coord_AVsh = x_sAV*cosd(gamma_stath)-y_sAV*sind(gamma_stath);

% MID 

cl_sm=teta_mid_stat/25;

% equations

y_sm = c_s*y_c_perc*cl_sm;
dy_dx_sm = dy_dx_perc*cl_sm;

eps_sm = atand(dy_dx_sm);

xss_sm = x_s-y_t_s.*sind(eps_sm);
yss_sm = y_sm+y_t_s.*cosd(eps_sm);

xps_sm = x_s+y_t_s.*sind(eps_sm);
yps_sm = y_sm-y_t_s.*cosd(eps_sm);

% Passage from x-y coordinates to T-ax

T_coord_ss_sm = -(xss_sm*sind(gamma_statm)+yss_sm*cosd(gamma_statm));
AX_coord_ss_sm = xss_sm*cosd(gamma_statm)-yss_sm*sind(gamma_statm);

T_coord_ps_sm = -(xps_sm*sind(gamma_statm)+yps_sm*cosd(gamma_statm));
AX_coord_ps_sm = xps_sm*cosd(gamma_statm)-yps_sm*sind(gamma_statm); 

T_coord_AVsm = -(x_sAV*sind(gamma_statm)+y_sAV*cosd(gamma_statm));
AX_coord_AVsm = x_sAV*cosd(gamma_statm)-y_sAV*sind(gamma_statm);

% TIP 

cl_st=teta_tip_stat/25;

% equations

y_st = c_s*y_c_perc*cl_st;
dy_dx_st = dy_dx_perc*cl_st;

eps_st = atand(dy_dx_st);

xss_st = x_s-y_t_s.*sind(eps_st);
yss_st = y_st+y_t_s.*cosd(eps_st);

xps_st = x_s+y_t_s.*sind(eps_st);
yps_st = y_st-y_t_s.*cosd(eps_st);

% Passage from x-y coordinates to T-ax

T_coord_ss_st = -(xss_st*sind(gamma_statt)+yss_st*cosd(gamma_statt));
AX_coord_ss_st = xss_st*cosd(gamma_statt)-yss_st*sind(gamma_statt);

T_coord_ps_st = -(xps_st*sind(gamma_statt)+yps_st*cosd(gamma_statt));
AX_coord_ps_st = xps_st*cosd(gamma_statt)-yps_st*sind(gamma_statt); 

T_coord_AVst = -(x_sAV*sind(gamma_statt)+y_sAV*cosd(gamma_statt));
AX_coord_AVst = x_sAV*cosd(gamma_statt)-y_sAV*sind(gamma_statt);

% HUB
DeltaT_AVsh = T_coord_AVsm - T_coord_AVsh;
DeltaAX_AVsh = AX_coord_AVsm - AX_coord_AVsh;

T_coord_ss_sh_trans = -(xss_sh*sind(gamma_stath)+yss_sh*cosd(gamma_stath))+DeltaT_AVsh;
AX_coord_ss_sh_trans = xss_sh*cosd(gamma_stath)-yss_sh*sind(gamma_stath)+DeltaAX_AVsh;

T_coord_ps_sh_trans = -(xps_sh*sind(gamma_stath)+yps_sh*cosd(gamma_stath))+DeltaT_AVsh;
AX_coord_ps_sh_trans = xps_sh*cosd(gamma_stath)-yps_sh*sind(gamma_stath)+DeltaAX_AVsh;

% TIP
DeltaT_AVst = T_coord_AVsm - T_coord_AVst;
DeltaAX_AVst = AX_coord_AVsm - AX_coord_AVst;

T_coord_ss_st_trans = -(xss_st*sind(gamma_statt)+yss_st*cosd(gamma_statt))+DeltaT_AVst;
AX_coord_ss_st_trans = xss_st*cosd(gamma_statt)-yss_st*sind(gamma_statt)+DeltaAX_AVst;

T_coord_ps_st_trans = -(xps_st*sind(gamma_statt)+yps_st*cosd(gamma_statt))+DeltaT_AVst;
AX_coord_ps_st_trans = xps_st*cosd(gamma_statt)-yps_st*sind(gamma_statt)+DeltaAX_AVst;

figure
hubsss = plot(T_coord_ss_sh_trans,-AX_coord_ss_sh_trans,'green');
hold on
hubsps = plot(T_coord_ps_sh_trans,-AX_coord_ps_sh_trans,'green');
midsss = plot(T_coord_ss_sm,-AX_coord_ss_rm,'red');
midsps = plot(T_coord_ps_sm,-AX_coord_ps_rm,'red');
tipsss = plot(T_coord_ss_st_trans,-AX_coord_ss_st_trans,'black');
tipsps = plot(T_coord_ps_st_trans,-AX_coord_ps_st_trans,'black');
title('STATOR BLADE (NACA 65)')
legend([hubsss midsss tipsss],{'HUB','MID','TIP'})
hold off

axis([-2.25 0.5 -2.5 0.25])

pbaspect([1 1 1])