

% Axial velocity @ rotor outlet (midspan) initialized to the inlet value

V_3A_m = linspace(V_2A_high(1)*0.8,V_2A_high(1),1000);

% Incidence calculation

alpha_2_geo  =  [alpha_2_geo_low(end:-1:1) alpha_2_geo_high(2:end) ];
i_2m = alpha_2_m - alpha_2_m_geo;
       
    for k = 1:length(V_3A_m)
    
        
    %% MID
    
    W_3A_m = V_3A_m(k);
    ddelta_di_S_m = ( 1+(sigma_S_m+0.25*sigma_S_m^4)*(abs(alpha_2_m)/53)^2.5 )/exp(3.1*sigma_S_m);
    delta_stat_od_m = dev_opt_statm + ddelta_di_S_m * (i_2m - i_opt_statm ) + 10*(1-V_3A_m(k)/V_2A_m);
    alpha_3_m = alpha_3_m_geo + delta_stat_od_m;
    V_3T_m = V_3A_m(k) * tand(alpha_3_m);
    V_3_m = sqrt(V_3T_m^2 + V_3A_m(k)^2);
    T_T3_m = T_T2_m;
    T_3_m = T_T2_m - V_3_m^2 / 2 / cp;
    
    Vmax_V2_m = 1.12 + 0.61 * (cosd(alpha_2_m))^2 / sigma_S_m * (V_2T_m - V_3T_m)/V_2A_m;  
    DmS = Vmax_V2_m * V_2_m / V_3_m;
    
    Y_3_p_tot_m_min = 0.004*(1+3.1*(DmS - 1)^2 + 0.4*(DmS - 1)^8)*2*sigma_S_m/cosd(alpha_3_m)* (V_3_m / V_2_m)^2;
    
    attack_S_c_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_c_m(end) - attack_S_c_m(end-1))/attack_S_c_m(end-1) ) > tol
        attack_S_c_m(end-1) = attack_S_c_m(end);
        
        alpha_c_m = attack_S_c_m(end) + gamma_statm;
        attack_S_c_m(end+1) = attack_S_m_design - 9 + ( 1 - (30/abs(alpha_c_m))^0.48 )*teta_mid_stat/4.176;
        
    end
    attack_S_c_m = attack_S_c_m(end);
    
    R_cm = attack_S_m_design - attack_S_c_m;
    
    i_c_m = i_opt_statm - R_cm / (1+0.5*M_2_m^3);
    
    % i_c_m = i_opt_statm + (attack_S_c_m - attack_S_m_design);
    
   attack_S_s_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_s_m(end) - attack_S_s_m(end-1))/attack_S_s_m(end-1) ) > tol
        attack_S_s_m(end-1) = attack_S_s_m(end);
        
        alpha_s_m = attack_S_s_m(end) + gamma_statm;
        attack_S_s_m(end+1) = attack_S_m_design + 10.3 + ( 2.92 - (abs(alpha_s_m)/15.6) )*teta_mid_stat/8.2;
        
    end
    attack_S_s_m = attack_S_s_m(end);
    
    R_sm = - attack_S_m_design + attack_S_s_m;
    
    i_s_m = i_opt_statm + R_sm/(1+0.5*(Ksh_i*M_2_m)^3);
    
    % i_s_m = i_opt_statm + (attack_S_s_m - attack_S_m_design);
    
    i_mm = i_c_m + (i_s_m - i_c_m)*R_cm/(R_cm + R_sm);
    
    % i_mm = i_min(alpha_2_m,alpha_3_m,th_c,sigma_S_m);
    
 %%% Critical Mach calculation
        
        T_cr_2_m = T_T2_m * 2 / (gamma+1);
        
        V_cr_2_m   = 1 * sqrt(gamma*R_star*T_cr_2_m);
        
        M_cr_2_m   = M_2_m * V_cr_2_m / ( V_2_m * Vmax_V2_m );
        
                if M_2_m > M_cr_2_m
        Y_3_p_tot_m_min = Y_3_p_tot_m_min + Ksh_i * ((M_2_m/M_cr_2_m-1)*V_cr_2_m/V_2_m)^2;
                end
    
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
    
 
    %% STATOR OUTLET high
    
    % Profiles
    V_3T_high(1) =  V_3T_m;
    p_3_high(1) = p_3_m;
    V_3A_high(1)    =  V_3A_m(k);
    V_3T_high(1)    =  V_3T_m(end);
    V_3_high(1)     =  V_3_m;
    alpha_3_high(1) = alpha_3_m;
    p_3_high(1)     = p_3_m;
    p_T3_high(1)    = p_T3_m;
    T_3_high(1)     = T_3_m;
    T_T3_high(1)    = T_T3_m;
    rho_3_high(1)   = rho_3_m;
    DS_high(1)      = DmS;

    i_2_high = alpha_2_high - alpha_2_geo_high;
    
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
    alpha_3_iter =  delta_stat_od_iter + alpha_3_geo_high(i);
    V_3A_iter(end+1) = V_3_iter * cosd(alpha_3_iter);
    
    end
    V_3A_iter = V_3A_iter(end);   
        
    V_3T_iter = V_3A_iter * tand(alpha_3_iter);
    V_3_iter = sqrt(V_3T_iter^2 + V_3A_iter(end)^2);
    
    Vmax_V2_iter = 1.12 + 0.61 * cosd(alpha_2_high(i))^2/sigma_S_high(i) * ( V_2T_high(i) - V_3T_iter ) / V_2A_high(i);
    Diter = Vmax_V2_iter * V_2_high(i) / V_3_iter;  
    
    Y_3_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_S_high(i) / cosd(alpha_3_iter) * (V_3_iter/V_2_high(i))^2;
    
    attack_S_c_iter = [attack_S_c_m*2 attack_S_c_m];
    while abs( (attack_S_c_iter(end) - attack_S_c_iter(end-1))/attack_S_c_iter(end-1) ) > tol
        attack_S_c_iter(end-1) = attack_S_c_iter(end);
        
        alpha_c_iter = attack_S_c_iter(end) + gamma_stat_high(i);
        
        if abs(alpha_c_iter) < 20
            
            alpha_c_iter = 20;
            
        end
       
        attack_S_c_iter(end+1) = attack_S_design_high(i) - 9 + ( 1 - (30/abs(alpha_c_iter))^0.48 )*teta_S_high(i)/4.176;
        
    end
    attack_S_c_iter = attack_S_c_iter(end);
    
    R_c_iter = attack_S_design_high(i) - attack_S_c_iter;
    
    i_c_iter = i_opt_stat_high(i) - R_c_iter / (1+0.5*M_2_high(i)^3);
    
    % i_c_iter = i_opt_stat_high(i) + (attack_S_c_iter - attack_S_design_high(i));
    
   attack_S_s_iter = [attack_S_s_m*2 attack_S_s_m];
    while abs( (attack_S_s_iter(end) - attack_S_s_iter(end-1))/attack_S_s_iter(end-1) ) > tol
        attack_S_s_iter(end-1) = attack_S_s_iter(end);
        
        alpha_s_iter = attack_S_s_iter(end) + gamma_stat_high(i);
        attack_S_s_iter(end+1) = attack_S_design_high(i) + 10.3 + ( 2.92 - (abs(alpha_s_iter)/15.6) )*teta_S_high(i)/8.2;
        
    end
    attack_S_s_iter = attack_S_s_iter(end);
    
    % i_s_iter = i_opt_stat_high(i) + (attack_S_s_iter - attack_S_design_high(i));
    
    R_s_iter = - attack_S_design_high(i) + attack_S_s_iter;
    
    i_s_iter = i_opt_stat_high(i) + R_s_iter/(1+0.5*(Ksh_i*M_2_high(i))^3);
    
    i_miter = i_c_iter + (i_s_iter - i_c_iter)*R_c_iter/(R_c_iter + R_s_iter);
    
    % i_miter = i_min(alpha_2_high(i),alpha_3_iter,th_c,sigma_S_high(i));
    
    %%% Critical Mach calculation
        
        T_cr_2_iter = T_T2_high(i) * 2 / (gamma+1);
        
        V_cr_2_iter   = 1 * sqrt(gamma*R_star*T_cr_2_iter);
      
        M_cr_2_iter   = M_2_high(i) * V_cr_2_iter / ( V_2_high(i) * Vmax_V2_iter );
        
                if M_R1_high(i) > M_cr_2_iter
        Y_3_p_tot_min_iter = Y_3_p_tot_min_iter + Ksh_i * ((M_2_high(i)/M_cr_2_iter-1)*V_cr_2_m/V_2_high(i))^2;
                end
    
    if i_2_high(i) >= i_miter
        
        xi_iter = (i_2_high(i) - i_miter) / (i_s_iter - i_miter);
        
    else
        
        xi_iter = (i_2_high(i) - i_miter) / (i_miter - i_c_iter);
        
    end
    
    if xi_iter <= 1 && xi_iter >= -2
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (1+xi_iter^2);
        
    elseif xi_iter < -2
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (5-4*(xi_iter + 2));
        
    else
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (2+2*(xi_iter-1));
        
    end
    
    p_T3_iter = p_T2_high(i) - Y_3_p_tot_iter * (p_T2_high(i) - p_2_high(i));
    p_3_iter  = p_3_high(i-1) / ( 1 -  V_3T_iter^2 / rhigh(i) / R_star / T_3_iter(end) * Drhigh );
    rho_3_iter = p_3_iter / R_star / T_3_iter(end);
    V_3_iter = sqrt( 2*gamma*R_star*T_3_iter(end) / (gamma-1) * ( (p_T3_iter/p_3_iter)^((gamma-1)/gamma) - 1 ) );
            
            T_3_iter(end+1) = T_T3_iter - V_3_iter^2/2/cp;
       
    end

    V_3A_high(i)    =  V_3A_iter;
    V_3T_high(i)    =  V_3T_iter;
    V_3_high(i)     =  V_3_iter;
    alpha_3_high(i) = alpha_3_iter;
    p_3_high(i)     = p_3_iter;
    p_T3_high(i)     = p_T3_iter;
    T_3_high(i)     = T_3_iter(end);
    rho_3_high(i)   = rho_3_iter;
    T_T3_high(i)   = T_T3_iter;
    DS_high(i)      = Diter;
    
    end
    
    %% STATOR OUTLET low
    
    % Profiles
    delta_od_IGV_low = zeros(1,length(rlow));
    V_3T_low(1) =  V_3T_m;
    p_3_low(1) = p_3_m;
    V_3A_low(1)    =  V_3A_m(k);
    V_3T_low(1)    =  V_3T_m(end);
    V_3_low(1)     =  V_3_m;
    alpha_3_low(1) = alpha_3_m;
    p_3_low(1)     = p_3_m;
    p_T3_low(1)    = p_T3_m;
    T_3_low(1)     = T_3_m;
    T_T3_low(1)   = T_T3_m;
    rho_3_low(1)   = rho_3_m;
    DS_low(1)      = DmS;

    i_2_low =  alpha_2_low - alpha_2_geo_low;
    
    s_S_low = pi * 2 * rlow / N_S;
    
    sigma_S_low = c_S_design ./ s_S_low;
    
    for i = 2:length(alpha_2_geo_low)
    
        T_3_iter = [T_3_low(i-1)*2 T_3_low(i-1)];
    while abs( (T_3_iter(end) - T_3_iter(end-1))/T_3_iter(end-1) ) > tol
        T_3_iter(end-1) = T_3_iter(end);
        
    T_T3_iter = T_T2_low(i);
    V_3_iter = sqrt(2*cp*(T_T3_iter-T_3_iter(end)));
    ddelta_di_S_iter = ( 1+(sigma_S_low(i)+0.25*sigma_S_low(i)^4)*(abs(alpha_2_low(i))/53)^2.5 )/exp(3.1*sigma_S_low(i));
    
        V_3A_iter = [V_3A_low(i-1)*2 V_3A_low(i-1)];
    while abs( (V_3A_iter(end) - V_3A_iter(end-1))/V_3A_iter(end-1) ) > tol
        V_3A_iter(end-1) = V_3A_iter(end);

    delta_stat_od_iter = dev_opt_stat_low(i) + ddelta_di_S_iter * (i_2_low(i) - i_opt_stat_low(i) ) + 10*(1-V_3A_iter(end)/V_2A_low(i));
    alpha_3_iter = delta_stat_od_iter + alpha_3_geo_low(i);
    V_3A_iter(end+1) = V_3_iter * cosd(alpha_3_iter);
    
    end
    V_3A_iter = V_3A_iter(end);   
        
    V_3T_iter = V_3A_iter * tand(alpha_3_iter);
    V_3_iter = sqrt(V_3T_iter^2 + V_3A_iter(end)^2);
    
    Vmax_V2_iter = 1.12 + 0.61 * cosd(alpha_2_low(i))^2/sigma_S_low(i) * ( V_2T_low(i) - V_3T_iter ) / V_2A_low(i);
    Diter = Vmax_V2_iter * V_2_low(i) / V_3_iter;  
    
    Y_3_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_S_low(i) / cosd(alpha_3_iter) * (V_3_iter/V_2_low(i))^2;
    
    attack_S_c_iter = [attack_S_c_m*2 attack_S_c_m];
    while abs( (attack_S_c_iter(end) - attack_S_c_iter(end-1))/attack_S_c_iter(end-1) ) > tol
        attack_S_c_iter(end-1) = attack_S_c_iter(end);
        
        alpha_c_iter = attack_S_c_iter(end) + gamma_stat_low(i);
        
        if abs(alpha_c_iter) < 20
            
            alpha_c_iter = 20;
            
        end
       
        attack_S_c_iter(end+1) = attack_S_design_low(i) - 9 + ( 1 - (30/abs(alpha_c_iter))^0.48 )*teta_S_low(i)/4.176;
        
    end
    attack_S_c_iter = attack_S_c_iter(end);
    
    R_c_iter = attack_S_design_low(i) - attack_S_c_iter;
    
    i_c_iter = i_opt_stat_low(i) - R_c_iter / (1+0.5*M_2_low(i)^3);
    
    % i_c_iter = i_opt_stat_low(i) + (attack_S_c_iter - attack_S_design_low(i));
    
   attack_S_s_iter = [attack_S_s_m*2 attack_S_s_m];
    while abs( (attack_S_s_iter(end) - attack_S_s_iter(end-1))/attack_S_s_iter(end-1) ) > tol
        attack_S_s_iter(end-1) = attack_S_s_iter(end);
        
        alpha_s_iter = attack_S_s_iter(end) + gamma_stat_low(i);
        attack_S_s_iter(end+1) = attack_S_design_low(i) + 10.3 + ( 2.92 - (abs(alpha_s_iter)/15.6) )*teta_S_low(i)/8.2;
        
    end
    attack_S_s_iter = attack_S_s_iter(end);
    
    % i_s_iter = i_opt_stat_low(i) + (attack_S_s_iter - attack_S_design_low(i));
    
    % i_miter = i_min(alpha_2_low(i),alpha_3_iter,th_c,sigma_S_low(i));
    
    R_s_iter = - attack_S_design_low(i) + attack_S_s_iter;
    
    i_s_iter = i_opt_stat_low(i) + R_s_iter/(1+0.5*(Ksh_i*M_2_low(i))^3);
    
    i_miter = i_c_iter + (i_s_iter - i_c_iter)*R_c_iter/(R_c_iter + R_s_iter);
    
    %%% Critical Mach calculation
        
        T_cr_2_iter = T_T2_low(i) * 2 / (gamma+1);
        
        V_cr_2_iter   = 1 * sqrt(gamma*R_star*T_cr_2_iter);
      
        M_cr_2_iter   = M_2_low(i) * V_cr_2_iter / ( V_2_low(i) * Vmax_V2_iter );
        
                if M_R1_low(i) > M_cr_2_iter
        Y_3_p_tot_min_iter = Y_3_p_tot_min_iter + Ksh_i * ((M_2_low(i)/M_cr_2_iter-1)*V_cr_2_m/V_2_low(i))^2;
                end
    
    if i_2_low(i) >= i_miter
        
        xi_iter = (i_2_low(i) - i_miter) / (i_s_iter - i_miter);
        
    else
        
        xi_iter = (i_2_low(i) - i_miter) / (i_miter - i_c_iter);
        
    end
    
    if xi_iter <= 1 && xi_iter >= -2
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (1+xi_iter^2);
        
    elseif xi_iter < -2
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (5-4*(xi_iter + 2));
        
    else
        
        Y_3_p_tot_iter = Y_3_p_tot_min_iter * (2+2*(xi_iter-1));
        
    end
    
    p_T3_iter = p_T2_low(i) - Y_3_p_tot_iter * (p_T2_low(i) - p_2_low(i));
    p_3_iter  = p_3_low(i-1) / ( 1 -  V_3T_iter^2 / rlow(i) / R_star / T_3_iter(end) * Drlow );
    rho_3_iter = p_3_iter / R_star / T_3_iter(end);
    V_3_iter = sqrt( 2*gamma*R_star*T_3_iter(end) / (gamma-1) * ( (p_T3_iter/p_3_iter)^((gamma-1)/gamma) - 1 ) );
            
            T_3_iter(end+1) = T_T3_iter - V_3_iter^2/2/cp;
       
    end

    V_3A_low(i)    =  V_3A_iter;
    V_3T_low(i)    =  V_3T_iter;
    V_3_low(i)     =  V_3_iter;
    alpha_3_low(i) = alpha_3_iter;
    p_3_low(i)     = p_3_iter;
    p_T3_low(i)     = p_T3_iter;
    T_3_low(i)     = T_3_iter(end);
    rho_3_low(i)   = rho_3_iter;
    T_T3_low(i)   = T_T3_iter;
    DS_low(i)      = Diter;
    
    end
    
    V_3A    =  [V_3A_low(end:-1:1) V_3A_high(2:end) ];
    V_3T    =  [V_3T_low(end:-1:1) V_3T_high(2:end) ];
    V_3     =  [V_3_low(end:-1:1) V_3_high(2:end) ];
    alpha_3 =  [alpha_3_low(end:-1:1) alpha_3_high(2:end) ];
    p_3     =  [p_3_low(end:-1:1) p_3_high(2:end) ];
    p_T3     =  [p_T3_low(end:-1:1) p_T3_high(2:end) ];
    T_3     =  [T_3_low(end:-1:1) T_3_high(2:end) ];
    T_T3   =  [T_T3_low(end:-1:1) T_T3_high(2:end) ];
    rho_3   =  [rho_3_low(end:-1:1) rho_3_high(2:end) ];
    i_2   =  [i_2_low(end:-1:1) i_2_high(2:end) ];
    DS      =  [DS_low(end:-1:1) DS_high(2:end)];

    
      dA = 2 * pi * (r(1:end-1)+r(2:end))/2 .* Dr;
      rho_3_av = (rho_3(1:end-1)+rho_3(2:end))/2;
      V_3A_av = (V_3A(1:end-1)+V_3A(2:end))/2;
      p_T3_ave = (p_T3(1:end-1)+p_T3(2:end))/2;
      T_T3_ave = (T_T3(1:end-1)+T_T3(2:end))/2;
      dm = rho_3_av .* dA .* V_3A_av;
      mnew(k) = sum(dm);
      
      if (mnew(k) > m*0.999) && (mnew(k) < m*1.001)
          V_3A_m = V_3A_m(k);
          break
      end
     
    end