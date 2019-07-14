        rho_2 = [rho_2_h(end) rho_2_m(end) rho_2_t(end)];
        V_2   = [V_2_h V_2_m V_2_t];
        W_2   = [W_2_h W_2_m W_2_t];
        T_2   = [T_2_h T_2_m T_2_t];
        p_2   = [p_2_h p_2_m p_2_t];
        
        MR_2   = W_2./sqrt(gamma*R_star*T_2);
        M_2    = V_2./sqrt(gamma*R_star*T_2);
          
        eta_R_av = mean([eta_R_h eta_R_m eta_R_t]);
        Y_p_R_av = mean([Y_2_p_tot_h Y_2_p_tot_m Y_2_p_tot_t]);
  