

% Axial velocity @ rotor outlet (midspan) initialized to the inlet value

V_3A_m = linspace(V_2A_high(1)*0.8,V_2A_high(1)*1.2,1000);

% Incidence calculation

alpha_2_geo  =  [alpha_2_geo_low(end:-1:1) alpha_2_geo_high(2:end) ];
i_2m = alpha_2_m_geo - alpha_2_m;
       
    for k = 1:length(V_3A_m)
    
        
    %% MID
    
    W_3A_m = V_3A_m(k);
    ddelta_di_S_m = ( 1+(sigma_S_m+0.25*sigma_S_m^4)*(abs(alpha_2_m)/53)^2.5 )/exp(3.1*sigma_S_m);
    delta_stat_od_m = dev_opt_statm + ddelta_di_S_m * (i_2m - i_opt_statm ) + 10*(1-V_3A_m(k)/V_2A_m);
    alpha_3_m = alpha_3_m_geo - delta_stat_od_m;
    V_3T_m = V_3A_m(k) * tand(alpha_3_m);
    W_3T_m = V_3T_m - U_m;
    W_3_m = sqrt(W_3T_m^2 + W_3A_m^2);
    V_3_m = sqrt(V_3T_m^2 + V_3A_m(k)^2);
    beta_3_m = atand(W_3T_m/W_3A_m(k));
    T_T3_m = T_T2_m;
    T_3_m = T_T2_m - V_3_m^2 / 2 / cp;
    
    Vmax_V2_m = 1.12 + 0.61 * (sind(alpha_2_m))^2 / sigma_S_m * (V_2T_m - V_3T_m)/V_2A_m;  
    DmS = Vmax_V2_m * V_2_m / V_3_m;
    
    Y_3_p_tot_m_min = 0.004*(1+3.1*(DmS - 1)^2 + 0.4*(DmS - 1)^8)*2*sigma_S_m/cosd(alpha_3_m)* (V_3_m / V_2_m)^2;
    
    attack_S_c_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_c_m(end) - attack_S_c_m(end-1))/attack_S_c_m(end-1) ) > tol
        attack_S_c_m(end-1) = attack_S_c_m(end);
        
        alpha_c_m = attack_S_c_m(end) + gamma_statm;
        attack_S_c_m(end+1) = attack_S_m_design - 9 + ( 1 - (30/abs(alpha_c_m))^0.48 )*teta_mid/4.176;
        
    end
    attack_S_c_m = attack_S_c_m(end);
    
    i_c_m = i_opt_rotm - (attack_S_c_m - attack_S_m_design);
    
   attack_S_s_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_s_m(end) - attack_S_s_m(end-1))/attack_S_s_m(end-1) ) > tol
        attack_S_s_m(end-1) = attack_S_s_m(end);
        
        alpha_s_m = attack_S_s_m(end) + gamma_statm;
        attack_S_s_m(end+1) = attack_S_m_design + 10.3 + ( 2.92 - (abs(alpha_s_m)/15.6) )*teta_mid/8.2;
        
    end
    attack_S_s_m = attack_S_s_m(end);
    
    i_s_m = i_opt_rotm + (attack_S_s_m - attack_S_m_design);
    
    i_mm = i_min(alpha_2_m,alpha_3_m,th_c,sigma_S_m);
    
    if i_2m >= i_mm
        
        xi_m = (i_2m - i_mm) / (i_s_m - i_mm);
        
    else
        
        xi_m = (i_2m - i_mm) / (i_mm - i_c_m);
        
    end
    
    if xi_m <= 1 && xi_m >= -2
        
        Y_3_p_tot_m = Y_3_p_tot_m_min * (1+xi_m^2);
        
    elseif xi_m < -2
        
        Y_3_p_tot_m = Y_3_p_tot_m_min * (5-4*(xi_m + 2));
        
    else
        
        Y_3_p_tot_m = Y_3_p_tot_m_min * (2+2*(xi_m-1));
        
    end
    
    p_T3_m = p_T2_m - Y_3_p_tot_m * (p_T2_m - p_2_m);
    p_3_m = p_T3_m / ( (1+(gamma-1)/2*(V_3_m^2/(gamma*R_star*T_3_m)))^(gamma/(gamma-1)) );
    rho_3_m = p_3_m / R_star / T_3_m;
    
 
    %% ROTOR INLET high
    
    % Profiles
    delta_od_IGV_high = zeros(1,length(rhigh));
    V_3T_high(1) =  V_3T_m;
    p_3_high(1) = p_3_m;
    V_3A_high(1)    =  V_3A_m(k);
    V_3T_high(1)    =  V_3T_m(end);
    V_3_high(1)     =  V_3_m;
    alpha_3_high(1) = alpha_3_m;
    W_3A_high(1)    =  W_3A_m;
    W_3T_high(1)    =  W_3T_m;
    W_3_high(1)     =  W_3_m;
    beta_3_high(1)  = beta_3_m;
    p_3_high(1)     = p_3_m;
    T_3_high(1)     = T_3_m;
    T_T3_high(1)   = T_T3_m;
    rho_3_high(1)   = rho_3_m;

    i_2_high = alpha_2_geo_high - alpha_2_high;
    
    s_S_high = pi * 2 * rhigh / N_S;
    
    sigma_S_high = c_S_design ./ s_S_high;
    
    for i = 2:length(alpha_2_geo_high)
    
        T_3_iter = [T_3_high(i-1)*2 T_3_high(i-1)];
    while abs( (T_3_iter(end) - T_3_iter(end-1))/T_3_iter(end-1) ) > tol
        T_3_iter(end-1) = T_3_iter(end);
        
    T_T3_iter = T_T2_high(i);
    V_3_iter = sqrt(2*cp*(T_T3_iter-T_3_iter(end)));
    ddelta_di_S_iter = ( 1+(sigma_S_high(i)+0.25*sigma_S_high(i)^4)*(abs(alpha_2_high(i))/53)^2.5 )/exp(3.1*sigma_S_high(i));
    
        V_3A_iter = [V_3A_high(i-1)*2 V_3A_high(i-1)];
    while abs( (V_3A_iter(end) - V_3A_iter(end-1))/V_3A_iter(end-1) ) > tol
        V_3A_iter(end-1) = V_3A_iter(end);

    delta_stat_od_iter = dev_opt_stat_high(i) + ddelta_di_S_iter * (i_2_high(i) - i_opt_stat_high(i) ) + 10*(1-V_3A_iter(end)/V_2A_high(i));
    beta_2_iter = beta_2_geo_high(i) - delta_rot_od_iter;
    V_3A_iter(end+1) = W_2_iter * cosd(beta_2_iter);
    
    end
    V_3A_iter = V_3A_iter(end);   
        
    V_2A_iter = W_2A_iter(end);
    W_2T_iter = W_2A_iter * tand(beta_2_iter);
    V_2T_iter = W_2T_iter + omega * rhigh(i);
    W_2_iter = sqrt(W_2T_iter^2 + W_2A_iter^2);
    V_2_iter = sqrt(V_2T_iter^2 + V_2A_iter(end)^2);
    alpha_2_iter = atand(V_2T_iter/V_2A_iter(end));
    
    Wmax_W1_iter = 1.12 + 0.61 * cosd(beta_1_high(i))^2/sigma_S_high(i) * ( V_1T_high(i) - V_2T_iter ) / V_1A_high(i);
    Diter = Wmax_W1_iter * W_1_high(i) / W_2_iter;  
    
    Y_2_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_S_high(i) / cosd(beta_2_iter) * (W_2_iter/W_1_high(i))^2;
    
    attack_c_iter = [attack_design_high(i)*2 attack_design_high(i)];
    while abs( (attack_c_iter(end) - attack_c_iter(end-1))/attack_c_iter(end-1) ) > tol
        attack_c_iter(end-1) = attack_c_iter(end);
        
        beta_c_iter = attack_c_iter(end) + gamma_rot_high(i);
        
        if abs(beta_c_iter) < 20
            
            beta_c_iter = 20;
            
        end
       
        attack_c_iter(end+1) = attack_design_high(i) - 9 + ( 1 - (30/abs(beta_c_iter))^0.48 )*teta_high(i)/4.176;
        
    end
    attack_c_iter = attack_c_iter(end);
    
    i_c_iter = i_opt_rot_high(i) - (attack_c_iter - attack_design_high(i));
    
   attack_s_iter = [attack_design_high(i)*2 attack_design_high(i)];
    while abs( (attack_s_iter(end) - attack_s_iter(end-1))/attack_s_iter(end-1) ) > tol
        attack_s_iter(end-1) = attack_s_iter(end);
        
        beta_s_iter = attack_s_iter(end) + gamma_rot_high(i);
        attack_s_iter(end+1) = attack_design_high(i) + 10.3 + ( 2.92 - (abs(beta_s_iter)/15.6) )*teta_high(i)/8.2;
        
    end
    attack_s_iter = attack_s_iter(end);
    
    i_s_iter = i_opt_rot_high(i) + (attack_s_iter - attack_design_high(i));
    
    i_miter = i_min(beta_1_high(i),beta_2_iter,th_c,sigma_S_high(i));
    
    if i_1_high(i) >= i_miter
        
        xi_iter = (i_1_high(i) - i_miter) / (i_s_iter - i_miter);
        
    else
        
        xi_iter = (i_1_high(i) - i_miter) / (i_miter - i_c_iter);
        
    end
    
    if xi_iter <= 1 && xi_iter >= -2
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (1+xi_iter^2);
        
    elseif xi_iter < -2
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (5-4*(xi_iter + 2));
        
    else
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (2+2*(xi_iter-1));
        
    end
    
    p_TR1_iter = p_1_high(i) * ( (1+(gamma-1)/2*(W_1_high(i)^2/(gamma*R_star*T_1_high(i))))^(gamma/(gamma-1)) );
    p_TR2_iter = p_TR1_iter - Y_2_p_tot_iter * (p_TR1_iter - p_1_high(i));
    p_2_iter  = p_2_high(i-1) / ( 1 -  V_2T_iter^2 / rhigh(i) / R_star / T_3_iter(end) * Drhigh );
    rho_2_iter = p_2_iter / R_star / T_3_iter(end);
    W_2_iter = sqrt( 2*gamma*R_star*T_3_iter(end) / (gamma-1) * ( (p_TR2_iter/p_2_iter)^((gamma-1)/gamma) - 1 ) );
            
            T_3_iter(end+1) = T_TR2_iter - W_2_iter^2/2/cp;
       
    end

    V_2A_high(i)    =  V_2A_iter;
    V_2T_high(i)    =  V_2T_iter;
    V_2_high(i)     =  V_2_iter;
    alpha_2_high(i) = alpha_2_iter;
    W_2A_high(i)    =  W_2A_iter;
    W_2T_high(i)    =  W_2T_iter;
    W_2_high(i)     =  W_2_iter;
    beta_2_high(i)  = beta_2_iter;
    p_2_high(i)     = p_2_iter;
    T_2_high(i)     = T_3_iter(end);
    rho_2_high(i)   = rho_2_iter;
    T_TR2_high(i)   = T_TR2_iter;
    
    end
    
    %% ROTOR INLET low
    
    % Profiles
    delta_od_IGV_low = zeros(1,length(rlow));
    V_2T_low(1)      =  V_2T_m;
    p_2_low(1)       = p_2_m;
    V_2A_low(1)      =  V_2A_m(k);
    V_2T_low(1)      =  V_2T_m(end);
    V_2_low(1)       =  V_2_m;
    alpha_2_low(1)   = alpha_2_m;
    W_2A_low(1)      =  W_2A_m;
    W_2T_low(1)      =  W_2T_m;
    W_2_low(1)       =  W_2_m;
    beta_2_low(1)    = beta_2_m;
    p_2_low(1)       = p_2_m;
    T_2_low(1)       = T_2_m;
    T_TR2_low(1)     = T_TR2_m;
    rho_2_low(1)     = rho_2_m;

    i_1_low = beta_1_geo_low - beta_1_low;
    
    s_R_low = pi * 2 * rlow / N_R;
    
    sigma_R_low = c_R_design ./ s_R_low;
    
    for i = 2:length(beta_1_geo_low)
    
        T_3_iter = [T_2_low(i-1)*2 T_2_low(i-1)];
    while abs( (T_3_iter(end) - T_3_iter(end-1))/T_3_iter(end-1) ) > tol
        T_3_iter(end-1) = T_3_iter(end);
        
    T_TR2_iter = T_TR1_low(i);
    W_2_iter = sqrt(2*cp*(T_TR2_iter-T_3_iter(end)));
    ddelta_di_iter = ( 1+(sigma_R_low(i)+0.25*sigma_R_low(i)^4)*(abs(beta_1_low(i))/53)^2.5 )/exp(3.1*sigma_R_low(i));
    
        W_2A_iter = [W_2A_low(i-1)*2 W_2A_low(i-1)];
    while abs( (W_2A_iter(end) - W_2A_iter(end-1))/W_2A_iter(end-1) ) > tol
        W_2A_iter(end-1) = W_2A_iter(end);

    delta_rot_od_iter = dev_opt_rot_low(i) + ddelta_di_iter * (i_1_low(i) - i_opt_rot_low(i) ) + 10*(1-W_2A_iter(end)/W_1A_low(i));
    beta_2_iter = beta_2_geo_low(i) - delta_rot_od_iter;
    W_2A_iter(end+1) = W_2_iter * cosd(beta_2_iter);
    
    end
    W_2A_iter = W_2A_iter(end);   
        
    V_2A_iter = W_2A_iter(end);
    W_2T_iter = W_2A_iter * tand(beta_2_iter);
    V_2T_iter = W_2T_iter + omega * rlow(i);
    W_2_iter = sqrt(W_2T_iter^2 + W_2A_iter^2);
    V_2_iter = sqrt(V_2T_iter^2 + V_2A_iter(end)^2);
    alpha_2_iter = atand(V_2T_iter/V_2A_iter(end));
    
    Wmax_W1_iter = 1.12 + 0.61 * cosd(beta_1_low(i))^2/sigma_R_low(i) * ( V_1T_low(i) - V_2T_iter ) / V_1A_low(i);
    Diter = Wmax_W1_iter * W_1_low(i) / W_2_iter;  
    
    Y_2_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_R_low(i) / cosd(beta_2_iter) * (W_2_iter/W_1_low(i))^2;
    
    attack_c_iter = [attack_design_low(i)*2 attack_design_low(i)];
    while abs( (attack_c_iter(end) - attack_c_iter(end-1))/attack_c_iter(end-1) ) > tol
        attack_c_iter(end-1) = attack_c_iter(end);
        
        beta_c_iter = attack_c_iter(end) + gamma_rot_low(i);
        
        if abs(beta_c_iter) < 20
            
            beta_c_iter = 20;
            
        end
        
        attack_c_iter(end+1) = attack_design_low(i) - 9 + ( 1 - (30/abs(beta_c_iter))^0.48 )*teta_low(i)/4.176;
        
    end
    attack_c_iter = attack_c_iter(end);
    
    i_c_iter = i_opt_rot_low(i) - (attack_c_iter - attack_design_low(i));
    
   attack_s_iter = [attack_design_low(i)*2 attack_design_low(i)];
    while abs( (attack_s_iter(end) - attack_s_iter(end-1))/attack_s_iter(end-1) ) > tol
        attack_s_iter(end-1) = attack_s_iter(end);
        
        beta_s_iter = attack_s_iter(end) + gamma_rot_low(i);
        attack_s_iter(end+1) = attack_design_low(i) + 10.3 + ( 2.92 - (abs(beta_s_iter)/15.6) )*teta_low(i)/8.2;
        
    end
    attack_s_iter = attack_s_iter(end);
    
    i_s_iter = i_opt_rot_low(i) + (attack_s_iter - attack_design_low(i));
    
    i_miter = i_min(beta_1_low(i),beta_2_iter,th_c,sigma_R_low(i));
    
    if i_1_low(i) >= i_miter
        
        xi_iter = (i_1_low(i) - i_miter) / (i_s_iter - i_miter);
        
    else
        
        xi_iter = (i_1_low(i) - i_miter) / (i_miter - i_c_iter);
        
    end
    
    if xi_iter <= 1 && xi_iter >= -2
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (1+xi_iter^2);
        
    elseif xi_iter < -2
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (5-4*(xi_iter + 2));
        
    else
        
        Y_2_p_tot_iter = Y_2_p_tot_min_iter * (2+2*(xi_iter-1));
        
    end
    
    p_TR1_iter = p_1_low(i) * ( (1+(gamma-1)/2*(W_1_low(i)^2/(gamma*R_star*T_1_low(i))))^(gamma/(gamma-1)) );
    p_TR2_iter = p_TR1_iter - Y_2_p_tot_iter * (p_TR1_iter - p_1_low(i));
    p_2_iter  = p_2_low(i-1) / ( 1 -  V_2T_iter^2 / rlow(i) / R_star / T_3_iter(end) * Drlow );
    rho_2_iter = p_2_iter / R_star / T_3_iter(end);
    W_2_iter = sqrt( 2*gamma*R_star*T_3_iter(end) / (gamma-1) * ( (p_TR2_iter/p_2_iter)^((gamma-1)/gamma) - 1 ) );
            
            T_3_iter(end+1) = T_TR2_iter - W_2_iter^2/2/cp;
       
    end

    V_2A_low(i)    =  V_2A_iter;
    V_2T_low(i)    =  V_2T_iter;
    V_2_low(i)     =  V_2_iter;
    alpha_2_low(i) = alpha_2_iter;
    W_2A_low(i)    =  W_2A_iter;
    W_2T_low(i)    =  W_2T_iter;
    W_2_low(i)     =  W_2_iter;
    beta_2_low(i)  = beta_2_iter;
    p_2_low(i)     = p_2_iter;
    T_2_low(i)     = T_3_iter(end);
    rho_2_low(i)   = rho_2_iter;
    T_TR2_low(i)   = T_TR2_iter;
    
    end
    
    V_2A    =  [V_2A_low(end:-1:1) V_2A_high(2:end) ];
    V_2T    =  [V_2T_low(end:-1:1) V_2T_high(2:end) ];
    V_2     =  [V_2_low(end:-1:1) V_2_high(2:end) ];
    alpha_2 =  [alpha_2_low(end:-1:1) alpha_2_high(2:end) ];
    W_2A    =  [W_2A_low(end:-1:1) W_2A_high(2:end) ];
    W_2T    =  [W_2T_low(end:-1:1) W_2T_high(2:end) ];
    W_2     =  [W_2_low(end:-1:1) W_2_high(2:end) ];
    beta_2  =  [beta_2_low(end:-1:1) beta_2_high(2:end) ];
    p_2     =  [p_2_low(end:-1:1) p_2_high(2:end) ];
    T_2     =  [T_2_low(end:-1:1) T_2_high(2:end) ];
    T_TR2   =  [T_TR2_low(end:-1:1) T_TR2_high(2:end) ];
    rho_2   =  [rho_2_low(end:-1:1) rho_2_high(2:end) ];

    
      dA = 2 * pi * (r(1:end-1)+r(2:end))/2 .* Dr;
      rho_2_av = (rho_1(1:end-1)+rho_1(2:end))/2;
      V_2A_av = (V_2A(1:end-1)+V_2A(2:end))/2;
      dm = rho_2_av .* dA .* V_2A_av;
      mnew(k) = sum(dm);
      
      if (mnew(k) > m*0.999) && (mnew(k) < m*1.001) 
          break
      end
     
    end