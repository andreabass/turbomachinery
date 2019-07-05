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
        
        mass_IGV   = [ m- rho_0(end) * V_0A * pi * b * D_m, m- mean(rho_1) * V_1A * pi * b * D_m ];
        energy_IGV = [ T_T0 - (T_1_t + V_1_t^2/2/cp) , T_T0 - (T_1_m + V_1_m^2/2/cp), T_T0 - (T_1_h + V_1_h^2/2/cp)];
       
    