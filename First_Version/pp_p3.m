        rho_3 = [rho_3_h(end) rho_3_m(end) rho_3_t(end)];
        V_3   = [V_3_h V_3_m V_3_t];
        T_3   = [T_3_h T_3_m T_3_t];
        p_3   = [p_3_h p_3_m p_3_t];
        
        M_3 = V_3./sqrt(gamma*R_star*T_3);
        
        eta_S_av = mean([eta_S_h eta_S_m eta_S_t]);
        Y_p_S_av = mean([Y_3_p_tot_h Y_3_p_tot_m Y_3_p_tot_t]);