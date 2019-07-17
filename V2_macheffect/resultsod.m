p_T3_av     = sum(p_T3_ave.*dm)/m;
beta_TT     = p_T3_av/p_T0;
T_T3_av     = sum(T_T3_ave.*dm)/m;
deltaHis_TT = cp * T_T0 * (beta_TT ^ ((gamma - 1)/gamma)-1); % [J/kg]
eta_TT      =  deltaHis_TT / ( cp*(T_T3_av - T_T0) );

Yigv    = (p_T0 - p_T1) ./ (p_T1 - p_1);
Yigv_av = sum(Yigv)/length(Yigv);

Yrot    = (p_TR1 - p_TR2) ./ (p_TR1 - p_1);
Yrot_av = sum(Yrot)/length(Yrot);

Ystat    = (p_T2 - p_T3) ./ (p_T2 - p_2);
Ystat_av = sum(Ystat)/length(Ystat);

eta_IGV = (T_0_m-T_1)./(T_0_m-T_0_m.*(p_1./p_0_m).^((gamma-1)/gamma));
eta_IGV_av = sum(eta_IGV)/length(eta_IGV);
eta_R = (T_1-T_1.*(p_2./p_1).^((gamma-1)/gamma)) ./ (T_1-T_2);
eta_R_av = sum(eta_R)/length(eta_R);
eta_S = (T_2-T_2.*(p_3./p_2).^((gamma-1)/gamma)) ./ (T_2-T_3);
eta_S_av = sum(eta_S)/length(eta_S);





        disp('_________________________________________________________________________________')
        disp(' ')
        disp('                             RESULTS (OFF-DESIGN)')
        disp('_________________________________________________________________________________')
        disp(' ')

disp(' ')
disp(['Total-to-total stage efficiency         :     ', num2str(eta_TT*100),'%'])
disp(['Static-to-static row efficiencies       :     ', num2str(eta_IGV_av*100), '%   ', num2str(eta_R_av*100), '%  ', num2str(eta_S_av*100),'%'])
disp(['Pressure loss coefficients              :     ', num2str(Yigv_av*100), '%    ', num2str(Yrot_av*100), '%   ', num2str(Ystat_av*100),'%'])
disp(' ')
disp(['Rotor diffusion coefficients   (h/m/t)  :     ', num2str([D(1) Dm D(end)])])
disp(['Stator diffusion coefficients  (h/m/t)  :     ', num2str([DS(1) DmS DS(end)])])
disp(' ')
disp(['IGV outlet Mach numbers        (h/m/t)  :     ', num2str(M_1(1)), '    ', num2str(M_1_m), '    ', num2str(M_1(end))])
disp(['Rotor inlet rel. Mach numbers  (h/m/t)  :     ', num2str(M_R1(1)), '    ', num2str(M_R1_m), '     ', num2str(M_R1(end))])
disp(['Stator inlet Mach numbers      (h/m/t)  :     ', num2str(M_2(1)), '    ', num2str(M_2_m), '    ', num2str(M_2(end))])



