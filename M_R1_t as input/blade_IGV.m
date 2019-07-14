%% NACA A4K6 SERIES BLADE DATA
% (suggested by AUNGIER)

x_perc=[0 1.25 2.5 5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percentage x/c
y_t_perc=[0 0.771 1.057 1.462 2.01 2.386 2.656 2.954 2.971 2.723 2.301 1.87 1.438 1.007 0.791 0];% half thickness t/c
y_c_perc=[0 0.792 1.357 2.248 3.531 4.42 5.04 5.71 5.82 5.516 4.891 4.011 2.922 1.642 0.912 0];% camber line y/c
dy_dx_perc=[Inf 0.5034 0.41 0.3131 0.2110 0.1483 0.1023 0.0359 0.0116 0.0478 0.0761 0.099 0.1184 0.1387 0.155 -Inf];% derivative

% THICKNESS PROFILE
x_IGV = x_perc*c_IGV;
y_t_IGV = c_IGV*y_t_perc*0.8;

%% HUB 

% CAMBER LINE PROFILE
cl_IGVh=teta_hub_IGV/25;
y_IGVh = c_IGV*y_c_perc*cl_IGVh;
dy_dx_IGVh = dy_dx_perc*cl_IGVh;
eps_IGVh = atand(dy_dx_IGVh);

% SUCTION SIDE PROFILE (X-Y)
xss_IGVh = x_IGV-y_t_IGV.*sind(eps_IGVh);
yss_IGVh = y_IGVh+y_t_IGV.*cosd(eps_IGVh);

% PRESSURE SIDE PROFILE (X-Y)
xps_IGVh = x_IGV+y_t_IGV.*sind(eps_IGVh);
yps_IGVh = y_IGVh-y_t_IGV.*cosd(eps_IGVh);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_IGVh = -(xss_IGVh*sind(gamma_IGVh)+yss_IGVh*cosd(gamma_IGVh));
AX_coord_ss_IGVh = xss_IGVh*cosd(gamma_IGVh)-yss_IGVh*sind(gamma_IGVh);

T_coord_ps_IGVh = -(xps_IGVh*sind(gamma_IGVh)+yps_IGVh*cosd(gamma_IGVh));
AX_coord_ps_IGVh = xps_IGVh*cosd(gamma_IGVh)-yps_IGVh*sind(gamma_IGVh);

%% MID 

% CAMBER LINE PROFILE
cl_IGVm=teta_mid_IGV/25;
y_IGVm = c_IGV*y_c_perc*cl_IGVm;
dy_dx_IGVm = dy_dx_perc*cl_IGVm;

% SUCTION SIDE PROFILE (X-Y)
eps_IGVm = atand(dy_dx_IGVm);
xss_IGVm = x_IGV-y_t_IGV.*sind(eps_IGVm);
yss_IGVm = y_IGVm+y_t_IGV.*cosd(eps_IGVm);

% PRESSURE SIDE PROFILE (X-Y)
xps_IGVm = x_IGV+y_t_IGV.*sind(eps_IGVm);
yps_IGVm = y_IGVm-y_t_IGV.*cosd(eps_IGVm);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_IGVm = -(xss_IGVm*sind(gamma_IGVm)+yss_IGVm*cosd(gamma_IGVm));
AX_coord_ss_IGVm = xss_IGVm*cosd(gamma_IGVm)-yss_IGVm*sind(gamma_IGVm);

T_coord_ps_IGVm = -(xps_IGVm*sind(gamma_IGVm)+yps_IGVm*cosd(gamma_IGVm));
AX_coord_ps_IGVm = xps_IGVm*cosd(gamma_IGVm)-yps_IGVm*sind(gamma_IGVm);

%% TIP 

% CAMBER LINE PROFILE
cl_IGVt=teta_tip_IGV/25;
y_IGVt = c_IGV*y_c_perc*cl_IGVt;
dy_dx_IGVt = dy_dx_perc*cl_IGVt;
eps_IGVt = atand(dy_dx_IGVt);

% SUCTION SIDE PROFILE (X-Y)
xss_IGVt = x_IGV-y_t_IGV.*sind(eps_IGVt);
yss_IGVt = y_IGVt+y_t_IGV.*cosd(eps_IGVt);

% PRESSURE SIDE PROFILE (X-Y)
xps_IGVt = x_IGV+y_t_IGV.*sind(eps_IGVt);
yps_IGVt = y_IGVt-y_t_IGV.*cosd(eps_IGVt);

% COORDINATE ROTATION (X-Y -> T-A)
T_coord_ss_IGVt = -(xss_IGVt*sind(gamma_IGVt)+yss_IGVt*cosd(gamma_IGVt));
AX_coord_ss_IGVt = xss_IGVt*cosd(gamma_IGVt)-yss_IGVt*sind(gamma_IGVt);

T_coord_ps_IGVt = -(xps_IGVt*sind(gamma_IGVt)+yps_IGVt*cosd(gamma_IGVt));
AX_coord_ps_IGVt = xps_IGVt*cosd(gamma_IGVt)-yps_IGVt*sind(gamma_IGVt); 