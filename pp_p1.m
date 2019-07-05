        rho_1 = [rho_1_h(end) rho_1_m(end) rho_1_t(end)];
        V_1   = [V_1_h V_1_m V_1_t];
        T_1   = [T_1_h T_1_m T_1_t];
        p_1   = [p_1_h p_1_m p_1_t];
        V_1T  = [V_1T_h V_1T_m(end) V_1T_t];
    
        story_rho_b   = b;
        story_rho_0   = rho_0;
        story_rho_1_m = rho_1_m;
        story_rho_1_t = rho_1_t; 
        story_rho_0   = rho_0; 
        story_V_1T_m  = V_1T_m;

        b = b(end);
        rho_0   = rho_0(end);
        rho_1_m = rho_1_m(end);
        rho_1_t = rho_1_t(end);
        rho_0   = rho_0(end);
        V_1T_m  = V_1T_m(end);
        
        T_1_m_is  = T_0_m * ( p_1_m / p_0_m )^((gamma-1)/gamma);
        T_1_t_is  = T_0_t * ( p_1_t / p_0_t )^((gamma-1)/gamma);
        T_1_h_is  = T_0_h * ( p_1_h / p_0_h )^((gamma-1)/gamma);
        
        eta_IGV_m = ( T_0_m - T_1_m ) / ( T_0_m - T_1_m_is );
        eta_IGV_t = ( T_0_t - T_1_t ) / ( T_0_t - T_1_t_is );
        eta_IGV_h = ( T_0_h - T_1_h ) / ( T_0_h - T_1_h_is );
        
        eta_IGV_av = mean([eta_IGV_t eta_IGV_m eta_IGV_h]);