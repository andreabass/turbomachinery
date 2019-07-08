    %% ROTOR INLET %%
     %%%%%%%%%%%%%%%
      %%%%%%%%%%%%%
       %%%%%%%%%%%
        %%%%%%%%%
         %%%%%%%
          %%%%%
           %%%
            %
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % Solution of simplified 1D problem as initial value for V_1T_m
        % (initialization of 1D problem: Vavra assumption for reaction degree).
        
        p1_1D
        
        V_1T_m = [2*V_1T(end) V_1T(end)];
       
    while abs((V_1T_m(end)-V_1T_m(end-1))/V_1T_m(end-1)) > tol
       
    V_1T_m(end-1) = V_1T_m(end);    
        
    %%%% MID %%%%

    V_1_m = sqrt(V_1T_m(end)^2 + V_1A_m^2);
    
    alpha_1_m = atand(V_1T_m(end) / V_1A_m);

    W_1T_m = V_1T_m(end) - U_m;

    W_1A_m = V_1A_m;

    W_1_m = sqrt(W_1A_m^2 + W_1T_m^2);

    beta_1_m = atand(W_1T_m / W_1A_m);
             
        alpha_2prime_m = 90 - alpha_1_m;
   
        if alpha_2prime_m > 30
        s_over_c_min_m = 0.614 + alpha_2prime_m / 130;
        else
   s_over_c_min_m = 0.46 + alpha_2prime_m / 77;
        end

        if alpha_2prime_m > 27
        A_m = 0.025 + (27 - alpha_2prime_m) / 3085;
        else
        A_m = 0.025 + (27 - alpha_2prime_m) / 530;
        end

        Re_ref = 2e5;
        Y_p_1_in = A_m;
        alpha_av_t_01 = atand((tand(alpha_0_m)+tand(alpha_1_m))/2);
        cL = 2 * s_over_c_min_m * (abs(tand(alpha_1_m)-tand(alpha_0_m)))*cosd(alpha_av_t_01);
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % IGV inlet density as initial value for rho_1_m
        
        rho_1_m = [2*rho_0(end) rho_0(end)]; 
        
    while abs((rho_1_m(end)-rho_1_m(end-1))/rho_1_m(end-1)) > tol
        
        rho_1_m(end-1) = rho_1_m(end); 
        
    Re_1_m = rho_1_m(end) * V_1_m * c_IGV / mu;
    
        Y_p_1_Re = Y_p_1_in * (Re_ref / Re_1_m)^0.2;
        Y_1_sec = c_IGV / b *(0.0334 * cosd(alpha_1_m)/cosd(alpha_0_m)) * (cL / s_over_c_min_m)^2 * ((cosd(alpha_1_m))^2) / (cosd(alpha_av_t_01))^3;
    
    Y_1_p_tot = Y_p_1_Re + Y_1_sec;
    
    T_T1_m = T_T0_m;
    
    T_1_m = T_T1_m - (V_1_m^2) / (2*cp);
    
    % p_T1_m = p_T0_m - Y_1_p_tot * (p_T0_m - p_0_m);
    
    p_T1_m = [2*p_T0_m p_T0_m];
    while abs((p_T1_m(end)-p_T1_m(end-1))/p_T1_m(end-1))>tol
    p_T1_m(end-1) = p_T1_m(end);
    p_1_m = p_T1_m(end) / ( (1+(gamma-1)/2*(V_1_m^2/(gamma*R_star*T_1_m)))^(gamma/(gamma-1)) );
    p_T1_m(end+1) = ( p_T0_m + Y_1_p_tot * p_1_m ) / (Y_1_p_tot+1);      
    end
    p_T1_m = p_T1_m(end); 

    rho_1_m(end+1) = p_1_m / R_star / T_1_m;
    
    end
       
    s_IGV_m = s_over_c_min_m * c_IGV;
    
    N_IGV = ceil(pi * D_m / s_IGV_m);
        s_IGV_m = pi * D_m / N_IGV;
        
    %%%% ROTOR INLET TIP %%%%
    
    s_IGV_t = pi * D_t / N_IGV;
    
    s_over_c_t = s_IGV_t / c_IGV;
    
    U_t = omega / 2 * D_t;
    
    V_1T_t = V_1T_m(end) * D_m / D_t; 

    alpha_1_t = atand(V_1T_t / V_1A_t);

    V_1_t = sqrt(V_1A_t^2 + V_1T_t^2);

    W_1T_t = V_1T_t - U_t;
    
    W_1A_t = V_1A_t;

    W_1_t = sqrt(W_1A_t^2 + W_1T_t^2);

    beta_1_t = atand(W_1T_t / W_1A_t);

        alpha_2prime_t = 90 - alpha_1_t;
        
        C_t = 0.08*((alpha_2prime_t/30)^2-1);
        n_AM_t = 1+alpha_2prime_t/30;
    
        if alpha_2prime_t > 27
        A_t = 0.025 + (27 - alpha_2prime_t) / 3085;
        else
        A_t = 0.025 + (27 - alpha_2prime_t) / 530;
        end
        
        if alpha_2prime_t > 30
        s_over_c_min_t = 0.614 + alpha_2prime_t / 130;
        B_t = 0; % Not available on Aungier (keep Y_p_1_in_t saturated at A)
        X_AM_t = s_over_c_t - s_over_c_min_t;
        Y_p_1_in_t = A_t + B_t * (abs(X_AM_t))^n_AM_t;
        else
        s_over_c_min_t = 0.46 + alpha_2prime_t / 77;
        B_t = 0.1583-alpha_2prime_t/1640;
        X_AM_t = s_over_c_t - s_over_c_min_t;
        Y_p_1_in_t = A_t+B_t*X_AM_t^2+C_t*X_AM_t^3;
        end

        alpha_av_t_01_t = atand((tand(alpha_0_t)+tand(alpha_1_t))/2);
        cL_t = 2 * s_over_c_t * (abs(tand(alpha_1_t)-tand(alpha_0_t)))*cosd(alpha_av_t_01_t);
    
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % Current density @ midspan as initial value for rho_1_t
        
        rho_1_t = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_t(end)-rho_1_t(end-1))/rho_1_t(end-1)) > tol
    
        rho_1_t(end-1) = rho_1_t(end); 
        
    Re_1_t = V_1_t * rho_1_t(end) * c_IGV / mu;
    
        Y_p_1_Re_t = Y_p_1_in_t * (Re_ref / Re_1_t)^0.2;
        Y_1_sec_t = c_IGV / b *(0.0334 * cosd(alpha_1_t)/cosd(alpha_0_t)) * (cL_t / s_over_c_t)^2 * ((cosd(alpha_1_t))^2) / (cosd(alpha_av_t_01_t))^3;
    
    Y_1_p_tot_t = Y_p_1_Re_t + Y_1_sec_t;
    
    T_T1_t = T_T0_t;

    T_1_t = T_T1_t - (V_1_t^2) / (2*cp);
    
    % p_T1_t = p_T0_t - Y_1_p_tot_t * (p_T0_t - p_0_t);
    
    p_T1_t = [2*p_T0_t p_T0_t];
    while abs((p_T1_t(end)-p_T1_t(end-1))/p_T1_t(end-1))>tol
    p_T1_t(end-1) = p_T1_t(end);
    p_1_t = p_T1_t(end) / ( (1+(gamma-1)/2*(V_1_t^2/(gamma*R_star*T_1_t)))^(gamma/(gamma-1)) );
    p_T1_t(end+1) = ( p_T0_t + Y_1_p_tot_t * p_1_t ) / (Y_1_p_tot_t+1);      
    end
    p_T1_t = p_T1_t(end);
    
    rho_1_t(end+1) = p_1_t / R_star / T_1_t;
    
    end

    %%%% ROTOR INLET HUB %%%%
    
    s_IGV_h = pi * D_h(end) / N_IGV;
    
    s_over_c_h = s_IGV_h / c_IGV;
    
    U_h = omega / 2 * D_h(end);

    V_1T_h = V_1T_m(end) * D_m / D_h(end); 

    alpha_1_h = atand(V_1T_h / V_1A_h);

    V_1_h = sqrt(V_1A_h^2 + V_1T_h^2);
    
    W_1A_h = V_1A_h;
    
    W_1T_h = V_1T_h - U_h;

    W_1_h = sqrt(W_1A_h^2 + W_1T_h^2);

    beta_1_h = atand(W_1T_h / W_1A_h);

        alpha_2prime_h = 90 - alpha_1_h;
        
        C_h = 0.08*((alpha_2prime_h/30)^2-1);
        n_AM_h = 1+alpha_2prime_h/30;
        
        if alpha_2prime_h > 27
        A_h = 0.025 + (27 - alpha_2prime_h) / 3085;
        else
        A_h = 0.025 + (27 - alpha_2prime_h) / 530;
        end
        
        if alpha_2prime_h > 30
        s_over_c_min_h = 0.614 + alpha_2prime_h / 130;
        B_h = 0; % Not available on Aungier (keep Y_p_1_in_t saturated at A)
        X_AM_h = s_over_c_h - s_over_c_min_h;
        Y_p_1_in_h = A_h + B_h * (abs(X_AM_h))^n_AM_h;
        else
        s_over_c_min_h = 0.46 + alpha_2prime_h / 77;
        B_h = 0.1583-alpha_2prime_h/1640;
        X_AM_h = s_over_c_h - s_over_c_min_h;
        Y_p_1_in_h = A_h+B_h*X_AM_h^2+C_h*X_AM_h^3;
        end
        
            alpha_av_t_01_h = atand((tand(alpha_0_h)+tand(alpha_1_h))/2);
            cL_h = 2 * s_over_c_h * (abs(tand(alpha_1_h)-tand(alpha_0_h)))*cosd(alpha_av_t_01_h);
    
    rho_1_h = 3 * rho_0(end) * V_0A / V_1A - rho_1_m(end) - rho_1_t(end);
        
    Re_1_h = V_1_h * rho_1_h * c_IGV / mu;
    
        Y_p_1_Re_h = Y_p_1_in_h * (Re_ref / Re_1_h)^0.2;
        Y_1_sec_h = c_IGV / b *(0.0334 * cosd(alpha_1_h)/cosd(alpha_0_h)) * (cL / s_over_c_h)^2 * ((cosd(alpha_1_h))^2) / (cosd(alpha_av_t_01_h))^3;
        
    Y_1_p_tot_h = Y_p_1_Re_h + Y_1_sec_h;
    
    T_T1_h = T_T0_h;

    T_1_h = T_T1_h - (V_1_h^2) / (2*cp);
    
    p_T1_h = [2*p_T0_h p_T0_h];
    while abs((p_T1_h(end)-p_T1_h(end-1))/p_T1_h(end-1))>tol
    p_T1_h(end-1) = p_T1_h(end);
    p_1_h = p_T1_h(end) / ( (1+(gamma-1)/2*(V_1_h^2/(gamma*R_star*T_1_h)))^(gamma/(gamma-1)) );
    p_T1_h(end+1) = ( p_T0_h + Y_1_p_tot_h * p_1_h ) / (Y_1_p_tot_h+1);      
    end
    p_T1_h = p_T1_h(end);
 
            V_1_h = sqrt( 2 * gamma * p_1_h / rho_1_h / (gamma-1) * ( (p_T1_h/p_1_h)^((gamma-1)/gamma) - 1 ) );
            
            V_1T_h = sqrt( V_1_h^2-V_1A^2 );
    
            V_1T_m(end+1) = V_1T_h * D_h / D_m;
    
    end
    
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