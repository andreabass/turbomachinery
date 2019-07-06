
figure

%STATOR
hubsss = plot(-T_coord_ss_sh_trans,-AX_coord_ss_sh_trans,'k--','LineWidth',1.5);
hold on
hubsps = plot(-T_coord_ps_sh_trans,-AX_coord_ps_sh_trans,'k--','LineWidth',1.5);
midsss = plot(-T_coord_ss_sm,-AX_coord_ss_sm,'k-.','LineWidth',1.5);
midsps = plot(-T_coord_ps_sm,-AX_coord_ps_sm,'k-.','LineWidth',1.5);
tipsss = plot(-T_coord_ss_st_trans,-AX_coord_ss_st_trans,'k','LineWidth',1.5);
tipsps = plot(-T_coord_ps_st_trans,-AX_coord_ps_st_trans,'k','LineWidth',1.5);

%ROTOR
hubrss = plot(T_coord_ss_rh_trans - T_coord_AVrt - T_coord_AVst,-AX_coord_ss_rh_trans + 3.2,'k--','LineWidth',1.5);
hubrps = plot(T_coord_ps_rh_trans - T_coord_AVrt - T_coord_AVst,-AX_coord_ps_rh_trans + 3.2,'k--','LineWidth',1.5);
midrss = plot(T_coord_ss_rm - T_coord_AVrt - T_coord_AVst,-AX_coord_ss_rm + 3.2,'k-.','LineWidth',1.5);
midrps = plot(T_coord_ps_rm - T_coord_AVrt - T_coord_AVst,-AX_coord_ps_rm + 3.2,'k-.','LineWidth',1.5);
tiprss = plot(T_coord_ss_rt_trans - T_coord_AVrt - T_coord_AVst,-AX_coord_ss_rt_trans + 3.2,'k','LineWidth',1.5);
tiprps = plot(T_coord_ps_rt_trans - T_coord_AVrt - T_coord_AVst,-AX_coord_ps_rt_trans + 3.2,'k','LineWidth',1.5);

%IGV
hubigvss = plot(-T_coord_ss_IGVh + 2,-AX_coord_ss_IGVh + 8.3,'k--','LineWidth',1.5);
hubigvps = plot(-T_coord_ps_IGVh + 2,-AX_coord_ps_IGVh + 8.3,'k--','LineWidth',1.5);
midigvss = plot(-T_coord_ss_IGVm + 2,-AX_coord_ss_IGVm + 8.3,'k-.','LineWidth',1.5);
midigvps = plot(-T_coord_ps_IGVm + 2,-AX_coord_ps_IGVm + 8.3,'k-.','LineWidth',1.5);
tipigvss = plot(-T_coord_ss_IGVt + 2,-AX_coord_ss_IGVt + 8.3,'k','LineWidth',1.5);
tipigvps = plot(-T_coord_ps_IGVt + 2,-AX_coord_ps_IGVt + 8.3,'k','LineWidth',1.5);

title('BLADE PROFILES')
legend([hubsss midsss tipsss],{'HUB','MID','TIP'})
hold off

axis([ -1.5 3.25 -3 8.875])

pbaspect([1 2.5 1])