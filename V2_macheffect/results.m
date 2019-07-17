

    p_3_av = 0.3*p_3_h + 0.4*p_3_m + 0.3*p_3_t;
    eta_TS = ( cp*T_0_m*( (p_3_av/p_0_m)^((gamma-1)/gamma) - 1 ) - V_0_m^2/2 ) /l_Eu;
    




        disp('_________________________________________________________________________________')
        disp(' ')
        disp('                                   RESULTS')
        disp('_________________________________________________________________________________')
        disp(' ')

disp(' ')
disp(['Total-to-total stage efficiency         :     ', num2str(eta_TT(end)*100),'%'])
disp(['Total-to-static stage efficiency        :     ', num2str(eta_TS*100),'%'])
disp(['Static-to-static row efficiencies       :     ', num2str(eta_IGV_av*100), '%   ', num2str(eta_R_av*100), '%  ', num2str(eta_S_av*100),'%'])
disp(['Pressure loss coefficients              :     ', num2str(Y_p_IGV_av*100), '%    ', num2str(Y_p_R_av*100), '%   ', num2str(Y_p_S_av*100),'%'])
disp(' ')
disp(['Tip peripheral velocity        [m/s]    :     ', num2str(U_m)]);
disp(' ')
disp(['Rotor diffusion coefficients   (h/m/t)  :     ', num2str([Dh Dm Dt])])
disp(['Stator diffusion coefficients  (h/m/t)  :     ', num2str([DhS DmS DtS])])
disp(' ')
disp(['IGV outlet Mach numbers        (h/m/t)  :     ', num2str(M_1(1)), '    ', num2str(M_1(2)), '    ', num2str(M_1(3))])
disp(['Rotor inlet rel. Mach numbers  (h/m/t)  :     ', num2str(MR_1(1)), '    ', num2str(MR_1(2)), '     ', num2str(MR_1(3))])
disp(['Stator inlet Mach numbers      (h/m/t)  :     ', num2str(M_2(1)), '    ', num2str(M_2(2)), '    ', num2str(M_2(3))])
disp(' ')
disp(['Number of blades                        :     ', num2str(N_IGV), '          ', num2str(N_R), '          ', num2str(N_S)])
Xi = cp * (T_2-T_1) / l_Eu;
disp(['Reaction degree                (h/m/t)  :     ', num2str(Xi(1)), '    ', num2str(Xi(2)), '     ', num2str(Xi(3))])
%disp(['According to Howell Statistics      :     ', Howell_R, Howell_S])




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

% STATOR
tvdstat = figure;
subplot(3,1,3)
velt(V_2A_h(end),V_2T_h,0)
velt(V_3A(end),V_3T_h,0,'k--')
title('STATOR HUB')
subplot(3,1,2)
velt(V_2A_m(end),V_2T_m,0)
velt(V_3A(end),V_3T_m,0,'k--')
title('STATOR MID')
subplot(3,1,1)
velt(V_2A_t(end),V_2T_t,0)
velt(V_3A(end),V_3T_t,0,'k--')
title('STATOR TIP')

%% OTHER OUTPUTS

Y_design      = [Y_p_IGV_av Y_p_R_av Y_p_S_av];
D_design_rot  = [Dh Dm Dt];
D_design_stat = [DhS DmS DtS];
eta_design   = [ eta_IGV_av eta_R_av eta_S_av ];