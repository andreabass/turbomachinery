
figure

XXstart = -1.5;
XXend = 6.25;
YYstart = -6.8;
YYend = 8.875;

TRASX_S = 4.5;
TRASY_S = -2.3;

%STATOR
hubsss = plot(T_coord_ss_sh_trans + TRASX_S,-AX_coord_ss_sh_trans + TRASY_S,'r');
hold on
hubsps = plot(T_coord_ps_sh_trans + TRASX_S,-AX_coord_ps_sh_trans + TRASY_S,'r','LineWidth',1.5);
midsss = plot(T_coord_ss_sm + TRASX_S,-AX_coord_ss_sm + TRASY_S,'g');
midsps = plot(T_coord_ps_sm + TRASX_S,-AX_coord_ps_sm + TRASY_S,'g');
tipsss = plot(T_coord_ss_st_trans + TRASX_S,-AX_coord_ss_st_trans + TRASY_S,'b');
tipsps = plot(T_coord_ps_st_trans + TRASX_S,-AX_coord_ps_st_trans + TRASY_S,'b');

TRASX_R = 0.5;
TRASY_R = 2.8;

%ROTOR
hubrss = plot(T_coord_ss_rh_trans - T_coord_AVrt - T_coord_AVst + TRASX_R,     -AX_coord_ss_rh_trans + TRASY_R,'r');
hubrps = plot(T_coord_ps_rh_trans - T_coord_AVrt - T_coord_AVst + TRASX_R,     -AX_coord_ps_rh_trans + TRASY_R,'r');
midrss = plot(T_coord_ss_rm - T_coord_AVrt - T_coord_AVst + TRASX_R,           -AX_coord_ss_rm + TRASY_R,'g');
midrps = plot(T_coord_ps_rm - T_coord_AVrt - T_coord_AVst + TRASX_R,           -AX_coord_ps_rm + TRASY_R,'g');
tiprss = plot(T_coord_ss_rt_trans - T_coord_AVrt - T_coord_AVst + TRASX_R,     -AX_coord_ss_rt_trans + TRASY_R,'b');
tiprps = plot(T_coord_ps_rt_trans - T_coord_AVrt - T_coord_AVst + TRASX_R,     -AX_coord_ps_rt_trans + TRASY_R,'b');

TRASX_IGV = 1;
TRASY_IGV = 8;

%IGV
hubigvss = plot(T_coord_ss_IGVh + TRASX_IGV, -AX_coord_ss_IGVh + TRASY_IGV,'r');
hubigvps = plot(T_coord_ps_IGVh + TRASX_IGV, -AX_coord_ps_IGVh + TRASY_IGV,'r');
midigvss = plot(T_coord_ss_IGVm + TRASX_IGV, -AX_coord_ss_IGVm + TRASY_IGV,'g');
midigvps = plot(T_coord_ps_IGVm + TRASX_IGV, -AX_coord_ps_IGVm + TRASY_IGV,'g');
tipigvss = plot(T_coord_ss_IGVt + TRASX_IGV, -AX_coord_ss_IGVt + TRASY_IGV,'b');
tipigvps = plot(T_coord_ps_IGVt + TRASX_IGV, -AX_coord_ps_IGVt + TRASY_IGV,'b');

title('BLADE PROFILES')
legend([hubsss midsss tipsss],{'HUB','MID','TIP'})
hold off

axis([ XXstart XXend YYstart YYend])

pbaspect([1 (YYend-YYstart)/(XXend-XXstart) 1])