      
    rho_1 = rho_0(end) * V_0A / V_1A;
    
        V_1T = U_m * ( 1 + 0.5 - deltaHis_TT / eta_TT(end) / U_m^2 / 2 );
        V_1T = [2*V_1T V_1T];
        
    while abs((V_1T(end)-V_1T(end-1))/V_1T(end-1)) > tol
        
        V_1T(end-1) = V_1T(end);

    V_1 = sqrt(V_1T(end)^2 + V_1A^2);
    
    alpha_1 = atand(V_1T(end) / V_1A);
    
    T_T1 = T_T0;

    T_1 = T_T1 - (V_1^2) / (2*cp);
    
    p_1 = R_star * rho_1 * T_1;
            
    Y_p_1_in = y_AM_inc_min(alpha_1);
        
    Re_1 = rho_1 * V_1 * c_IGV / mu;
    
    Y_p_1_Re = y_AM_Re(Y_p_1_in,Re_1); 
    
    s_over_c_min = s_c_min_AM(alpha_1);

    Y_1_sec = y_AM_sec(alpha_1,alpha_0_m,c_IGV,b,s_over_c_min);
        
    Y_1_p_tot = Y_p_1_Re + Y_1_sec;
    
    p_T1 = ( p_T0 + Y_1_p_tot * p_1 ) / (Y_1_p_tot+1); 
    
    V_1 = sqrt( 2*gamma*R_star*T_1 / (gamma-1) * ( (p_T1/p_1)^((gamma-1)/gamma) - 1 ) );
       
    s_IGV = s_over_c_min * c_IGV;
    
    N_bl_IGV = ceil(pi * D_m / s_IGV);
    
        s_IGV = pi * D_m / N_bl_IGV;
    
        alpha_1 = acosd(V_1A/V_1);
    
        V_1T(end+1) = V_1 * sind(alpha_1);
    
    end
