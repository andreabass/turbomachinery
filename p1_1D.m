      
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
             
        alpha_2prime = 90 - alpha_1;
   
        if alpha_2prime > 30
        s_over_c_min = 0.614 + alpha_2prime / 130;
        else
   s_over_c_min = 0.46 + alpha_2prime / 77;
        end

        if alpha_2prime > 27
        A = 0.025 + (27 - alpha_2prime) / 3085;
        else
        A = 0.025 + (27 - alpha_2prime) / 530;
        end

        Re_ref = 2e5;
        Y_p_1_in = A;
        alpha_av_t_01 = atand((tand(alpha_0_m)+tand(alpha_1))/2);
        cL = 2 * s_over_c_min * (abs(tand(alpha_1)-tand(alpha_0_m)))*cosd(alpha_av_t_01);
   
    Re_1 = rho_1(end) * V_1 * c_IGV / mu;
   
        Y_p_1_Re = Y_p_1_in * (Re_ref / Re_1)^0.2;
        Y_1_sec = c_IGV / b *(0.0334 * cosd(alpha_1)/cosd(alpha_0_m)) * (cL / s_over_c_min)^2 * ((cosd(alpha_1))^2) / (cosd(alpha_av_t_01))^3;
    
    Y_1_p_tot = Y_p_1_Re + Y_1_sec;
    
    p_T1 = ( p_T0 + Y_1_p_tot * p_1 ) / (Y_1_p_tot+1); 
    
    V_1 = sqrt( 2*gamma*R_star*T_1 / (gamma-1) * ( (p_T1/p_1)^((gamma-1)/gamma) - 1 ) );
       
    s_IGV = s_over_c_min * c_IGV;
    
    N_bl_IGV = ceil(pi * D_m / s_IGV);
    
        s_IGV = pi * D_m / N_bl_IGV;
    
        alpha_1 = acosd(V_1A/V_1);
    
        V_1T(end+1) = V_1 * sind(alpha_1);
    
    end
