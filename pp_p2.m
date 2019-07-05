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
        
        mass_R   = [ m - mean(rho_1) * V_1A * pi * b *D_m, m - mean(rho_2) * V_2A * pi * b * D_m ];
        energy_R = [ l_Eu - cp*(T_T2_m - T_T1_m) ; l_Eu - cp*(T_T2_t - T_T1_t); l_Eu - cp*(T_T2_h - T_T1_h)];
       
    