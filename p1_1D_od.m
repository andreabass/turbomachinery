     

        %%% [ INITIALIZATION ]  %%%
        % rho_1 initialized to the design value of the density @ rotor mid
        % N.B. the IGV deviation angle is neglected (it is just an
        % initialization)
        
        rho_1 = [2*rho_1_m rho_1_m];
        
   while abs((rho_1(end)-rho_1(end-1))/rho_1(end-1)) > tol
        
   rho_1(end-1) = rho_1(end);
        
    V_1A = rho_0(end) * V_0A / rho_1(end);
    
    W_1 = V_1A / cosd(beta_1_m);
    
    W_1T = W_1 * sind(beta_1_m);
    
    V_1T = W_1T + U_m;

    V_1 = sqrt(V_1T^2 + V_1A^2);
    
    alpha_1 = atand(V_1T / V_1A);
    
    IGV_rotation = alpha_1 - alpha_1_m_geo; % IGV deviation angle is neglected
    
    i_IGV_od = IGV_rotation; % This because the inlet flow is axial and the design incidence angle is null
    
    T_T1 = T_T0;

    T_1 = T_T1 - (V_1^2) / (2*cp);
    
    p_1 = R_star * rho_1(end) * T_1;
    
    alpha_0_geo = i_IGV_od;
    
    Y_p_1_tot = y_AM_od(alpha_0_geo, alpha_0_m, alpha_1, sigma_IGV_m, i_IGV_od, rho_1(end), b, c_IGV, V_1, mu);
    
    p_T1 = ( p_T0 + Y_1_p_tot * p_1 ) / (Y_1_p_tot+1); 
    
    p_1 = p_T1 / ( 1 + V_1^2 / ( 2*gamma*R_star*T_1 / (gamma-1) ) )^( gamma/(gamma-1) );
        
        rho_1(end+1) = p_1 / R_star / T_1;
    
   end
