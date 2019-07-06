disp(' ')
disp(['Total-to-total stage efficiency        :     ', num2str([eta_TT(end)])])
disp(['Static-to-static row efficiencies      :     ', num2str([eta_IGV_av eta_R_av eta_S_av])])
disp(['Pressure loss coefficients             :     ', num2str([Y_p_IGV_av Y_p_R_av Y_p_S_av])])
disp(['Rotor diffusion coefficients (h/m/t)   :     ', num2str([Dh Dm Dt])])
disp(['Stator diffusion coefficients (h/m/t)  :     ', num2str([DhS DmS DtS])])
disp(['Reaction degree at midspan             :     ', num2str(Xi_m)])
%disp(['According to Howell Statistics      :     ', Howell_R, Howell_S])

