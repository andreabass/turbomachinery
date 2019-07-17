

% Axial velocity @ rotor outlet (midspan) initialized to the inlet value

V_2A_m = linspace(V_1A_high(1)*0.8,V_1A_high(1)*1.2,1000);

% Incidence calculation

beta_1_geo  =  [beta_1_geo_low(end:-1:1) beta_1_geo_high(2:end) ];
i_1m = beta_1_m_od - beta_1_m_geo;

% Mach number @ inlet

M_R1_m = W_1_m / sqrt(gamma*R_star*T_1_m);
       
    for k = 1:length(V_2A_m)
    
        
    %% MID
    
    W_2A_m = V_2A_m(k);
    ddelta_di_m = ( 1+(sigma_R_m+0.25*sigma_R_m^4)*(abs(beta_1_m_od)/53)^2.5 )/exp(3.1*sigma_R_m);
    delta_rot_od_m = dev_opt_rotm + ddelta_di_m * (i_1m - i_opt_rotm ) + 10*(1-W_2A_m/W_1A_m);
    beta_2_m = beta_2_m_geo - delta_rot_od_m;
    W_2T_m = W_2A_m * tand(beta_2_m);
    V_2T_m = W_2T_m + U_m;
    W_2_m = sqrt(W_2T_m^2 + W_2A_m^2);
    V_2_m = sqrt(V_2T_m^2 + V_2A_m(k)^2);
    alpha_2_m = atand(V_2T_m/V_2A_m(k));
    T_TR2_m = T_TR1_m;
    T_2_m = T_TR2_m - W_2_m^2 / 2 / cp;
    T_T2_m = T_2_m + V_2_m^2/2/cp;
    
    Km = abs(tand(beta_1_m) - V_2A(end)/V_1A_m * tand(beta_2_m));
    
    Wmax_W1_m = 1.12 + 0.61 * cosd(beta_1_m)^2/sigma_R_m * Km;
    Dm = Wmax_W1_m * W_1_m / W_2_m;  
    
    Y_2_p_tot_m_min = 0.004 * ( 1 + 3.1*(Dm-1)^2 + 0.4*(Dm-1)^8 ) * 2 * sigma_R_m / cosd(beta_2_m) * (W_2_m/W_1_m)^2;
    
    attack_c_m = [attack_m_design*2 attack_m_design];
    while abs( (attack_c_m(end) - attack_c_m(end-1))/attack_c_m(end-1) ) > tol
        attack_c_m(end-1) = attack_c_m(end);
        
        beta_c_m = -attack_c_m(end) + gamma_rotm;
        attack_c_m(end+1) = attack_m_design - 9 + ( 1 - (30/abs(beta_c_m))^0.48 )*teta_mid/4.176;
        
    end
    attack_c_m = attack_c_m(end);
    
    R_cm = attack_m_design - attack_c_m;
    
    i_c_m = i_opt_rotm - R_cm / (1+0.5*M_R1_m^3);
    
   attack_s_m = [attack_m_design*2 attack_m_design];
    while abs( (attack_s_m(end) - attack_s_m(end-1))/attack_s_m(end-1) ) > tol
        attack_s_m(end-1) = attack_s_m(end);
        
        beta_s_m = -attack_s_m(end) + gamma_rotm;
        attack_s_m(end+1) = attack_m_design + 10.3 + ( 2.92 - (abs(beta_s_m)/15.6) )*teta_mid/8.2;
        
    end
    attack_s_m = attack_s_m(end);
    
    R_sm = - attack_m_design + attack_s_m;
    
    i_s_m = i_opt_rotm + R_sm/(1+0.5*(Ksh_i*M_R1_m)^3);
    
    % i_mm = i_min(beta_1_m,beta_2_m,th_c,sigma_R_m);
    
    i_mm = i_c_m + (i_s_m - i_c_m)*R_cm/(R_cm + R_sm);
    
    % Y_2_p_tot_m_min = Y_2_p_tot_m_min * (1+(i_mm-i_opt_rotm)^2/R_sm);
    
    %%% Critical Mach calculation
    
        T_TR1_m = T_1_m + W_1_m^2/2/cp;
        
        T_cr_1_m = T_TR1_m * 2 / (gamma+1);
        
        W_cr_1_m   = 1 * sqrt(gamma*R_star*T_cr_1_m);
        
        M_R1_m     = W_1_m / sqrt(gamma*R_star*T_1_m);
        
        M_cr_1_m   = M_R1_m * W_cr_1_m / ( W_1_m * Wmax_W1_m );
        
                if M_R1_m > M_cr_1_m
        Y_2_p_tot_m_min = Y_2_p_tot_m_min + Ksh_i * ((M_R1_m/M_cr_1_m-1)*W_cr_1_m/W_1_m)^2;
                end
    
    
    if i_1m >= i_mm
        
        xi_m = (i_1m - i_mm) / (i_s_m - i_mm);
        
    else
        
        xi_m = (i_1m - i_mm) / (i_mm - i_c_m);
        
    end
    
    if xi_m <= 1 && xi_m >= -2
        
        Y_2_p_tot_m = Y_2_p_tot_m_min * (1+xi_m^2);
        
    elseif xi_m < -2
        
        Y_2_p_tot_m = Y_2_p_tot_m_min * (5-4*(xi_m + 2));
        
    else
        
        Y_2_p_tot_m = Y_2_p_tot_m_min * (2+2*(xi_m-1));
        
    end
    
    
    p_TR1_m = p_1_m * ( (1+(gamma-1)/2*(W_1_m^2/(gamma*R_star*T_1_m)))^(gamma/(gamma-1)) );
    p_TR2_m = p_TR1_m - Y_2_p_tot_m * (p_TR1_m - p_1_m);
    p_2_m = p_TR2_m / ( (1+(gamma-1)/2*(W_2_m^2/(gamma*R_star*T_2_m)))^(gamma/(gamma-1)) );
    rho_2_m = p_2_m / R_star / T_2_m;
    p_T2_m = p_2_m * ( (1+(gamma-1)/2*(rho_2_m*V_2_m^2/(gamma*p_2_m)))^(gamma/(gamma-1)) );
    
 
    %% ROTOR OUTLET high
    
    % Profiles
    delta_od_IGV_high = zeros(1,length(rhigh));
    V_2T_high(1) =  V_2T_m;
    p_2_high(1) = p_2_m;
    V_2A_high(1)    =  V_2A_m(k);
    V_2T_high(1)    =  V_2T_m(end);
    V_2_high(1)     =  V_2_m;
    alpha_2_high(1) = alpha_2_m;
    W_2A_high(1)    =  W_2A_m;
    W_2T_high(1)    =  W_2T_m;
    W_2_high(1)     =  W_2_m;
    beta_2_high(1)  = beta_2_m;
    p_2_high(1)     = p_2_m;
    T_2_high(1)     = T_2_m;
    T_TR2_high(1)   = T_TR2_m;
    rho_2_high(1)   = rho_2_m;
    p_T2_high(1)    = p_T2_m;
    p_TR1_high(1)   = p_TR1_m;
    p_TR2_high(1)   = p_TR2_m;
    D_high(1)       = Dm;

    i_1_high = beta_1_high - beta_1_geo_high;
    
    s_R_high = pi * 2 * rhigh / N_R;
    
    sigma_R_high = c_R_design ./ s_R_high;
    
    for i = 2:length(beta_1_geo_high)
    
        T_2_iter = [T_2_high(i-1)*2 T_2_high(i-1)];
    while abs( (T_2_iter(end) - T_2_iter(end-1))/T_2_iter(end-1) ) > tol
        T_2_iter(end-1) = T_2_iter(end);
        
    T_TR2_iter = T_TR1_high(i);
    W_2_iter = sqrt(2*cp*(T_TR2_iter-T_2_iter(end)));
    ddelta_di_iter = ( 1+(sigma_R_high(i)+0.25*sigma_R_high(i)^4)*(abs(beta_1_high(i))/53)^2.5 )/exp(3.1*sigma_R_high(i));
    
        W_2A_iter = [W_2A_high(i-1)*2 W_2A_high(i-1)];
    while abs( (W_2A_iter(end) - W_2A_iter(end-1))/W_2A_iter(end-1) ) > tol
        W_2A_iter(end-1) = W_2A_iter(end);

    delta_rot_od_iter = dev_opt_rot_high(i) + ddelta_di_iter * (i_1_high(i) - i_opt_rot_high(i) ) + 10*(1-W_2A_iter(end)/W_1A_high(i));
    beta_2_iter = beta_2_geo_high(i) - delta_rot_od_iter;
    W_2A_iter(end+1) = W_2_iter * cosd(beta_2_iter);
    
    end
    W_2A_iter = W_2A_iter(end);   
        
    V_2A_iter = W_2A_iter(end);
    W_2T_iter = W_2A_iter * tand(beta_2_iter);
    V_2T_iter = W_2T_iter + omega * rhigh(i);
    W_2_iter = sqrt(W_2T_iter^2 + W_2A_iter^2);
    V_2_iter = sqrt(V_2T_iter^2 + V_2A_iter(end)^2);
    alpha_2_iter = atand(V_2T_iter/V_2A_iter(end));
    
    Kiter = abs(tand(beta_1_high(i)) - V_2A_iter/V_1A_high(i) * tand(beta_2_iter));
    
    Wmax_W1_iter = 1.12 + 0.61 * cosd(beta_1_high(i))^2/sigma_R_high(i) * Kiter;
    Diter = Wmax_W1_iter * W_1_high(i) / W_2_iter;  
    
    Y_2_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_R_high(i) / cosd(beta_2_iter) * (W_2_iter/W_1_high(i))^2;
    
    attack_c_iter = [attack_design_high(i)*2 attack_design_high(i)];
    while abs( (attack_c_iter(end) - attack_c_iter(end-1))/attack_c_iter(end-1) ) > tol
        attack_c_iter(end-1) = attack_c_iter(end);
        
        beta_c_iter = -attack_c_iter(end) + gamma_rot_high(i);
        
        if abs(beta_c_iter) < 30

            beta_c_iter = 30;
            
        end
       
        attack_c_iter(end+1) = attack_design_high(i) - 9 + ( 1 - (30/abs(beta_c_iter))^0.48 )*teta_high(i)/4.176;
        
    end
    attack_c_iter = attack_c_iter(end);
    
    R_c_iter = attack_design_high(i) - attack_c_iter;
    
    i_c_iter = i_opt_rot_high(i) - R_c_iter / (1+0.5*M_R1_high(i)^3);
    
    % i_c_iter = i_opt_rot_high(i) + (attack_c_iter - attack_design_high(i));
    
   attack_s_iter = [attack_design_high(i)*2 attack_design_high(i)];
    while abs( (attack_s_iter(end) - attack_s_iter(end-1))/attack_s_iter(end-1) ) > tol
        attack_s_iter(end-1) = attack_s_iter(end);
        
        beta_s_iter = -attack_s_iter(end) + gamma_rot_high(i);
        attack_s_iter(end+1) = attack_design_high(i) + 10.3 + ( 2.92 - (abs(beta_s_iter)/15.6) )*teta_high(i)/8.2;
        
    end
    attack_s_iter = attack_s_iter(end);
    
    R_s_iter = - attack_design_high(i) + attack_s_iter;
    
    i_s_iter = i_opt_rot_high(i) + R_s_iter/(1+0.5*(Ksh_i*M_R1_high(i))^3);
    
    % i_s_iter = i_opt_rot_high(i) + (attack_s_iter - attack_design_high(i));
    
    % i_miter = i_min(beta_1_high(i),beta_2_iter,th_c,sigma_R_high(i));
    
    i_miter = i_c_iter + (i_s_iter - i_c_iter)*R_c_iter/(R_c_iter + R_s_iter);
    
    %%% Critical Mach calculation
        
        T_cr_1_iter = T_TR1_high(i) * 2 / (gamma+1);
        
        W_cr_1_iter   = 1 * sqrt(gamma*R_star*T_cr_1_iter);
      
        M_cr_1_iter   = M_R1_high(i) * W_cr_1_iter / ( W_1_high(i) * Wmax_W1_iter );
        
                if M_R1_high(i) > M_cr_1_iter
        Y_2_p_tot_min_iter = Y_2_p_tot_min_iter + Ksh_i * ((M_R1_high(i)/M_cr_1_iter-1)*W_cr_1_m/W_1_high(i))^2;
                end
    
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
    p_2_iter  = p_2_high(i-1) / ( 1 -  V_2T_iter^2 / rhigh(i) / R_star / T_2_iter(end) * Drhigh );
    rho_2_iter = p_2_iter / R_star / T_2_iter(end);
    W_2_iter = sqrt( 2*gamma*R_star*T_2_iter(end) / (gamma-1) * ( (p_TR2_iter/p_2_iter)^((gamma-1)/gamma) - 1 ) );
    p_T2_iter = p_2_iter * ( (1+(gamma-1)/2*(rho_2_iter*V_2_iter^2/(gamma*p_2_iter)))^(gamma/(gamma-1)) );
            
            T_2_iter(end+1) = T_TR2_iter - W_2_iter^2/2/cp;
       
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
    T_2_high(i)     = T_2_iter(end);
    rho_2_high(i)   = rho_2_iter;
    T_TR2_high(i)   = T_TR2_iter;    
    p_T2_high(i)    = p_T2_iter;
    p_TR1_high(i)    = p_TR1_iter;
    p_TR2_high(i)    = p_TR2_iter;
    D_high(i)       = Diter;
    
    end
    
    %% ROTOR OUTLET low
    
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
    p_T2_low(1)      = p_T2_m;
    p_TR1_low(1)    = p_TR1_m;
    p_TR2_low(1)    = p_TR2_m;
    K_low(1)        = Km;
    D_low(1)       = Dm;
    
    
    i_c_iter_low(1) = i_c_m;
    i_miter_low(1) = i_mm;
    beta_c_low(1) = beta_c_m;

    i_1_low =  beta_1_low - beta_1_geo_low;
    
    s_R_low = pi * 2 * rlow / N_R;
    
    sigma_R_low = c_R_design ./ s_R_low;
    
    for i = 2:length(beta_1_geo_low)
    
        T_2_iter = [T_2_low(i-1)*2 T_2_low(i-1)];
    while abs( (T_2_iter(end) - T_2_iter(end-1))/T_2_iter(end-1) ) > tol
        T_2_iter(end-1) = T_2_iter(end);
        
    T_TR2_iter = T_TR1_low(i);
    W_2_iter = sqrt(2*cp*(T_TR2_iter-T_2_iter(end)));
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

    Kiter = abs(tand(beta_1_low(i)) - V_2A_iter/V_1A_low(i) * tand(beta_2_iter));
    
    Wmax_W1_iter = 1.12 + 0.61 * cosd(beta_1_low(i))^2/sigma_R_low(i) * Kiter;
    
    Diter = Wmax_W1_iter * W_1_low(i) / W_2_iter;  
    
    Y_2_p_tot_min_iter = 0.004 * ( 1 + 3.1*(Diter-1)^2 + 0.4*(Diter-1)^8 ) * 2 * sigma_R_low(i) / cosd(beta_2_iter) * (W_2_iter/W_1_low(i))^2;
    
    attack_c_iter = [attack_design_low(i)*2 attack_design_low(i)];
    while abs( (attack_c_iter(end) - attack_c_iter(end-1))/attack_c_iter(end-1) ) > tol
        attack_c_iter(end-1) = attack_c_iter(end);
        
        beta_c_iter = -attack_c_iter(end) + gamma_rot_low(i);
        
        if abs(beta_c_iter) < 20
            
            beta_c_iter = 20;
            
        end
        
        attack_c_iter(end+1) = attack_design_low(i) - 9 + ( 1 - (30/abs(beta_c_iter))^0.48 )*teta_low(i)/4.176;
        
    end
    attack_c_iter = attack_c_iter(end);
    
    R_c_iter = attack_design_low(i) - attack_c_iter;
    
    i_c_iter = i_opt_rot_low(i) - R_c_iter / (1+0.5*M_R1_low(i)^3);
    
    % i_c_iter = i_opt_rot_low(i) + (attack_c_iter - attack_design_low(i));
    
   attack_s_iter = [attack_design_low(i)*2 attack_design_low(i)];
    while abs( (attack_s_iter(end) - attack_s_iter(end-1))/attack_s_iter(end-1) ) > tol
        attack_s_iter(end-1) = attack_s_iter(end);
        
        beta_s_iter = -attack_s_iter(end) + gamma_rot_low(i);
        attack_s_iter(end+1) = attack_design_low(i) + 10.3 + ( 2.92 - (abs(beta_s_iter)/15.6) )*teta_low(i)/8.2;
        
    end
    attack_s_iter = attack_s_iter(end);
    
    i_s_iter = i_opt_rot_low(i) + (attack_s_iter - attack_design_low(i));
    
    R_s_iter = - attack_design_low(i) + attack_s_iter;
    
    i_s_iter = i_opt_rot_low(i) + R_s_iter/(1+0.5*(Ksh_i*M_R1_low(i))^3);
    
    i_miter = i_c_iter + (i_s_iter - i_c_iter)*R_c_iter/(R_c_iter + R_s_iter);
    
    % i_miter = i_min(beta_1_low(i),beta_2_iter,th_c,sigma_R_low(i));
    
 %%% Critical Mach calculation
        
        T_cr_1_iter = T_TR1_low(i) * 2 / (gamma+1);
        
        W_cr_1_iter   = 1 * sqrt(gamma*R_star*T_cr_1_iter);
      
        M_cr_1_iter   = M_R1_low(i) * W_cr_1_iter / ( W_1_low(i) * Wmax_W1_iter );
        
                if M_R1_low(i) > M_cr_1_iter
        Y_2_p_tot_min_iter = Y_2_p_tot_min_iter + Ksh_i * ((M_R1_low(i)/M_cr_1_iter-1)*W_cr_1_m/W_1_low(i))^2;
                end

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
    p_2_iter  = p_2_low(i-1) / ( 1 -  V_2T_iter^2 / rlow(i) / R_star / T_2_iter(end) * Drlow );
    rho_2_iter = p_2_iter / R_star / T_2_iter(end);
    W_2_iter = sqrt( 2*gamma*R_star*T_2_iter(end) / (gamma-1) * ( (p_TR2_iter/p_2_iter)^((gamma-1)/gamma) - 1 ) );
    p_T2_iter = p_2_iter * ( (1+(gamma-1)/2*(rho_2_iter*V_2_iter^2/(gamma*p_2_iter)))^(gamma/(gamma-1)) );
            
        T_2_iter(end+1) = T_TR2_iter - W_2_iter^2/2/cp;
       
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
    T_2_low(i)     = T_2_iter(end);
    rho_2_low(i)   = rho_2_iter;
    T_TR2_low(i)   = T_TR2_iter;    
    p_T2_low(i)    = p_T2_iter;
    p_TR1_low(i)    = p_TR1_iter;
    p_TR2_low(i)    = p_TR2_iter;
    K_low(i)        = Kiter;
    D_low(i)       = Diter;
    
    i_c_iter_low(i) = i_c_iter;
    i_miter_low(i) = i_miter;
    beta_c_low(i) = beta_c_iter;
    
    end
    
    T_T2_high    = T_2_high + V_2_high.^2/2/cp;
    T_T2_low    = T_2_low + V_2_low.^2/2/cp;
    
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
    i_1     =  [i_1_low(end:-1:1) i_1_high(2:end)];
    T_T2    =  [T_T2_low(end:-1:1) T_T2_high(2:end)];
    p_T2    =  [p_T2_low(end:-1:1) p_T2_high(2:end)];
    p_TR2   =  [p_TR2_low(end:-1:1) p_TR2_high(2:end)];
    p_TR1   =  [p_TR1_low(end:-1:1) p_TR1_high(2:end)];
    D       =  [D_low(end:-1:1) D_high(2:end)];

    
      dA = 2 * pi * (r(1:end-1)+r(2:end))/2 .* Dr;
      rho_2_av = (rho_2(1:end-1)+rho_2(2:end))/2;
      V_2A_av = (V_2A(1:end-1)+V_2A(2:end))/2;
      dm = rho_2_av .* dA .* V_2A_av;
      mnew(k) = sum(dm);
      
      if (mnew(k) > m*0.999) && (mnew(k) < m*1.001)
          V_2A_m = V_2A_m(k);
          break
      end
   
    end
    
    M_2_high = V_2_high ./ sqrt(gamma*R_star.*T_2_high);
    M_2_low = V_2_low ./ sqrt(gamma*R_star.*T_2_low);
    M_2     =  [M_2_low(end:-1:1) M_2_high(2:end) ];
    M_2_m   = V_2_m / sqrt(gamma*R_star*T_2_m);
    
    
    
