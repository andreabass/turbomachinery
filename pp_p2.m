        rho_2 = [rho_2_h(end) rho_2_m(end) rho_2_t(end)];
        V_2   = [V_2_h V_2_m V_2_t];
        T_2   = [T_2_h T_2_m T_2_t];
        p_2   = [p_2_h p_2_m p_2_t];
    
        story_eta_R_t = eta_R_t;
        story_eta_R_m = eta_R_m; 
        story_eta_R_h = eta_R_h; 
        story_V_2A    = V_2A;

        eta_R_m   = eta_R_m(end);
        eta_R_t   = eta_R_t(end);
        eta_R_h   = eta_R_h(end);
        V_2A      = V_2A(end);
        
        eta_R_av = mean([eta_R_t eta_R_m eta_R_h]);
        Y_p_R_av = mean([Y_2_p_tot_h Y_2_p_tot_m Y_2_p_tot_t]);
  