
    V_0T_m = 0;
    V_0T_t = 0; 
    V_0T_h = 0;
    
    phi_1_m = 0.5;
          
    b =  (D_t - D_h)/2 ;
    D_m = (D_h + D_t)/2;
    U_m = omega * D_m / 2;
    
    V_1A_m     = phi_1_m * omega * D_m / 2;
          V_1A       = V_1A_m; % (definition)
    
    V_1A_t = V_1A_m;
    V_1A_h = V_1A_m;
  
    c_IGV = 0.04;
    
    c_IGV_t = c_IGV;
    c_IGV_h = c_IGV;
  
    t_over_s_h = 0.02;
    
    %% IGV INLET %%
    %%%%%%%%%%%%%%%
     %%%%%%%%%%%%%
      %%%%%%%%%%%
       %%%%%%%%%
        %%%%%%%
         %%%%%
          %%%
           %
    
        %%%%%% [INITIALIZATION] %%%%%%%
      
        % Total density @ IGV inlet as initial value for rho_0
        
        rho_0 = [rho_T0 rho_T0+2*tol ];

        while abs(rho_0(end)-rho_0(end-1)) > tol 
        
        rho_0(end-1) = rho_0(end);    
        
    V_0A = m / ( pi * b * D_m * rho_0(end) ); % (eq. 15)

        V_0A_m   = V_0A; % (definition)
    V_0A_t = V_0A_m;     % (eq. 16)
    V_0A_h = V_0A_m;     % (eq. 17)

    alpha_0_t = atand(V_0T_t / V_0A_t); % (eq. 18)
    alpha_0_m = atand(V_0T_m / V_0A_m); % (eq. 19)
    alpha_0_h = atand(V_0T_h / V_0A_h); % (eq. 20)

    V_0_t = sqrt(V_0A_t^2 + V_0T_t^2); % (eq. 21)
    V_0_m = sqrt(V_0A_m^2 + V_0T_m^2); % (eq. 22)
    V_0_h = sqrt(V_0A_h^2 + V_0T_h^2); % (eq. 23)

    T_0_t = T_T0_t - (V_0_t^2)/(2*cp); % (eq. 24)
    T_0_m = T_T0_m - (V_0_m^2)/(2*cp); % (eq. 25)
    T_0_h = T_T0_h - (V_0_h^2)/(2*cp); % (eq. 26)

    p_0_t = p_T0_t / (1 + (V_0_t^2)/(2 * R_star * T_0_t)); % (eq. 27)
    p_0_m = p_T0_m / (1 + (V_0_m^2)/(2 * R_star * T_0_m)); % (eq. 28)
    p_0_h = p_T0_h / (1 + (V_0_h^2)/(2 * R_star * T_0_h)); % (eq. 29)

    rho_0_t = p_0_t / R_star / T_0_t; % (eq. 30)
    rho_0_m = p_0_m / R_star / T_0_m; % (eq. 31)
    rho_0_h = p_0_h / R_star / T_0_h; % (eq. 32)
        rho_0(end+1) = rho_0_m;              % (definition)
        
        end

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
        % and rho_1_m (initialization of 1D problem: Vavra assumption for
        % reaction degree).
        
        p1_1D
        
        V_1T_m = [V_1T V_1T+2*tol];
        rho_1_m = [rho_1 rho_1+2*tol];
        b = (D_t - D_h) / 2;
        b = [b b+2*tol];
       
    while abs(rho_1_m(end)-rho_1_m(end-1)) > tol || abs(V_1T_m(end)-V_1T_m(end-1)) > tol
        
        rho_1_m(end-1) = rho_1_m(end);  
        V_1T_m(end-1) = V_1T_m(end);
        
        while abs(b(end)-b(end-1))>tol
        b(end-1)   = b(end);
        D_m        = D_t - b(end);
   V_1A_m     = phi_1_m * omega * D_m / 2;
        V_1A       = V_1A_m; % (definition)
        b(end+1)   = m / ( rho_1_m(end) * pi * D_m * V_1A);
        end
        
    D_h = D_t - 2*b(end);
    D_m = (D_h + D_t)/2;
    U_m = omega * D_m / 2;
        
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
        
        % Current density @ midspan as initial value for rho_1_t
        
        rho_1_m = [rho_1_m(end) rho_1_m(end)+2*tol]; 
        
    while abs(rho_1_m(end)-rho_1_m(end-1)) > tol
    
        rho_1_m(end-1) = rho_1_m(end); 
        
    Re_1_m = rho_1_m(end) * V_1_m * c_IGV / mu;
    
        Y_p_1_Re = Y_p_1_in * (Re_ref / Re_1_m)^0.2;
        Y_1_sec = c_IGV / b(end) *(0.0334 * cosd(alpha_1_m)/cosd(alpha_0_m)) * (cL / s_over_c_min_m)^2 * ((cosd(alpha_1_m))^2) / (cosd(alpha_av_t_01))^3;
    
    Y_1_p_tot = Y_p_1_Re + Y_1_sec;
    
    p_T1_m = p_T0_m - Y_1_p_tot * (p_T0_m - p_0_m);
    
    T_T1_m = T_T0_m;
    
    T_1_m = T_T1_m - (V_1_m^2) / (2*cp);
    
    p_1_m = p_T1_m / (1 + (V_1_m^2)/(2 * R_star * T_1_m));
    
    rho_1_m(end+1) = p_1_m / R_star / T_1_m;
    
    end
       
    s_IGV_m = s_over_c_min_m * c_IGV;
    
    N_bl_IGV = pi * D_m / s_IGV_m;
        
    %%%% ROTOR INLET TIP %%%%
    
    s_IGV_t = pi * D_t / N_bl_IGV;
    
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
        
        rho_1_t = [rho_1_m(end) rho_1_m(end)+2*tol]; 
    
    while abs(rho_1_t(end)-rho_1_t(end-1)) > tol
    
        rho_1_t(end-1) = rho_1_t(end); 
        
    Re_1_t = V_1_t * rho_1_t(end) * c_IGV / mu;
    
        Y_p_1_Re_t = Y_p_1_in_t * (Re_ref / Re_1_t)^0.2;
        Y_1_sec_t = c_IGV / b(end) *(0.0334 * cosd(alpha_1_t)/cosd(alpha_0_t)) * (cL_t / s_over_c_t)^2 * ((cosd(alpha_1_t))^2) / (cosd(alpha_av_t_01_t))^3;
    
    Y_1_p_tot_t = Y_p_1_Re_t + Y_1_sec_t;
    
    T_T1_t = T_T0_t;

    T_1_t = T_T1_t - (V_1_t^2) / (2*cp);
    
    p_T1_t = p_T0_t - Y_1_p_tot_t * (p_T0_t - p_0_t);
    
    p_1_t = p_T1_t / (1 + (V_1_t^2)/(2 * R_star * T_1_t));
    
    rho_1_t(end+1) = p_1_t / R_star / T_1_t;
    
    end

    %%%% ROTOR INLET HUB %%%%
    
    s_IGV_h = pi * D_h / N_bl_IGV;
    
    s_over_c_h = s_IGV_h / c_IGV;
    
    U_h = omega / 2 * D_h;

    V_1T_h = V_1T_m(end) * D_m / D_h; 

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
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % Current density @ midspan as initial value for rho_1_h
        
            rho_1_h = [rho_1_m(end) rho_1_m(end)+2*tol]; 
            alpha_av_t_01_h = atand((tand(alpha_0_h)+tand(alpha_1_h))/2);
            cL_h = 2 * s_over_c_h * (abs(tand(alpha_1_h)-tand(alpha_0_h)))*cosd(alpha_av_t_01_h);
    
    while abs(rho_1_h(end)-rho_1_h(end-1)) > tol
    
        rho_1_h(end-1) = rho_1_h(end); 
        
    Re_1_h = V_1_h * rho_1_h(end) * c_IGV / mu;
    
        Y_p_1_Re_h = Y_p_1_in_h * (Re_ref / Re_1_h)^0.2;
        Y_1_sec_h = c_IGV / b(end) *(0.0334 * cosd(alpha_1_h)/cosd(alpha_0_h)) * (cL / s_over_c_h)^2 * ((cosd(alpha_1_h))^2) / (cosd(alpha_av_t_01_h))^3;
        
    Y_1_p_tot_h = Y_p_1_Re_h + Y_1_sec_h;
    
    T_T1_h = T_T0_h;

    T_1_h = T_T1_h - (V_1_h^2) / (2*cp);
    
    p_T1_h = p_T0_h - Y_1_p_tot_h * (p_T0_h - p_0_h);
    
    p_1_h = p_T1_h / (1 + (V_1_h^2)/(2 * R_star * T_1_h));
    
    rho_1_h(end+1) = p_1_h / R_star / T_1_h;
    
    end
    
    rho_1_m(end+1) = 3 * rho_0(end) * V_0A / V_1A - rho_1_h(end) - rho_1_t(end);
    
            V_1_m = sqrt( 2*(p_T1_m - p_1_m)/rho_1_m(end) );
            
            alpha_1_m = acosd(V_1A/V_1_m);
    
            V_1T_m(end) = V_1_m * sind(alpha_1_m);
    
    end
 
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
        V_1T_m  = V_1T_m(end);
        
        %Mass and energy conservation check and initialization variables
        %convergence check
        mass   = [ rho_0(end) * V_0A * pi * b *D_m, mean(rho_1) * V_1A * pi * b * D_m ];
        energy = [ T_T0 T_1_t + V_1_t^2/2/cp ; T_T0 T_1_m + V_1_m^2/2/cp; T_T0 T_1_h + V_1_h^2/2/cp];
        
        figure(1)
        plot(story_rho_b);
        figure(2)
        plot(story_rho_0);
        figure(3)
        plot(story_rho_1_m);
        figure(4)
        plot(story_rho_1_t);
        figure(5)
        plot(story_V_1T_m);
        
        
        %% ROTOR  OUTLET %%
          %%%%%%%%%%%%%%%
           %%%%%%%%%%%%%
            %%%%%%%%%%%
             %%%%%%%%%
              %%%%%%%
               %%%%%
                %%%
                 %
        
                 
        %%%%%% [INITIALIZATION] %%%%%%%
        %Initialization value for V_2A_m
        V_2A_m = V_1A_m;
        V_2A_t = V_2A_m;
        V_2A_h = V_2A_m;
        
        % Let's define initialization values for inlet axial velocity and
        % rotor efficiency defined as the ratio between isentropic static
        % enthalpy change and real static enthalpy change. Thanks to these
        % values we'll be able to fully determine the VT and TDN quantities
        % at the rotor outlet that will allow us to calculate losses and
        % therefore get to adjusted values of the TDN quantities. Such
        % updated values will be finally used to recalculate the efficiency
        % of the rotor.In the meanwhile for every cycle we'll impose the
        % mass balance to update the axial velocity.
        
        V_2A = [V_2A_m V_2A_m + 2*tol];
        eta_R_t = [eta_TT_m - 2*tol eta_TT_m];
        eta_R_m = [eta_TT_m - 2*tol eta_TT_m];
        eta_R_h = [eta_TT_m - 2*tol eta_TT_m];
        
        while abs(eta_R_t(end) - eta_R_t(end-1))> tol || abs(eta_R_m(end) - eta_R_m(end-1))> tol || abs(eta_R_h(end) - eta_R_h(end-1))> tol
        while abs(V_2A(end) - V_2A(end-1))> tol
        V_2A(end-1) = V_2A;
        
        l_Eu = deltaHis_TT / eta_TT_m;
        V_2T_m = V_1T_m - l_Eu / U_m;
        V_2T_t = V_2T_m * D_m / D_t;
        V_2T_h = V_2T_m * D_m / D_h;
        
        W_2T_m = V_2T_m - U_m;
        W_2T_t = V_2T_t - U_t;
        W_2T_h = V_2T_h - U_h;
        
        W_2A_m = V_2A(end);
        W_2A_t = V_2A(end);
        W_2A_h = V_2A(end);
        
        W_2A = W_2A_m;
        
        V_2_t = sqrt(V_2T_t^2 + V_2A(end)^2);
        V_2_m = sqrt(V_2T_m^2 + V_2A(end)^2);
        V_2_h = sqrt(V_2T_h^2 + V_2A(end)^2);
        
        W_2_t = sqrt(W_2T_t^2 + W_2A^2);
        W_2_m = sqrt(W_2T_m^2 + W_2A^2);
        W_2_h = sqrt(W_2T_h^2 + W_2A^2);
        
        alpha_2_t = atand(V_2T_t / V_2A(end));
        alpha_2_m = atand(V_2T_m / V_2A(end));
        alpha_2_h = atand(V_2T_h / V_2A(end));
        
        beta_2_t = atand(W_2T_t / W_2A);
        beta_2_m = atand(W_2T_m / W_2A);
        beta_2_h = atand(W_2T_h / W_2A);
        
        Xi_t = (W_1_t^2 - W_2_t^2) / 2 / l_Eu;
        Xi_m = (W_1_m^2 - W_2_m^2) / 2 / l_Eu;
        Xi_h = (W_1_h^2 - W_2_h^2) / 2 / l_Eu;
        
        DeltaH_R_t = Xi_t * l_Eu;
        DeltaH_R_m = Xi_m * l_Eu;
        DeltaH_R_h = Xi_h * l_Eu;
        
        T_2_t = T_1_t + DeltaH_R_t / cp;
        T_2_m = T_1_m + DeltaH_R_m / cp;
        T_2_h = T_1_h + DeltaH_R_h / cp;
        
        T_2_t_is = eta_R_t(end) * (T_2_t - T_1_t);
        T_2_m_is = eta_R_m(end) * (T_2_m - T_1_m);
        T_2_h_is = eta_R_h(end) * (T_2_h - T_1_h);
        
        % The outlet pressure is equal to the one we'd reach with an
        % isentropic process
        
        p_2_t = p_1_t * (T_1_t / T_2_t_is) ^ (gamma/(1-gamma));
        p_2_m = p_1_m * (T_1_m / T_2_m_is) ^ (gamma/(1-gamma));
        p_2_h = p_1_h * (T_1_h / T_2_h_is) ^ (gamma/(1-gamma));
        
        rho_2_t = p_2_t / R_star / T_2_t;
        rho_2_m = p_2_m / R_star / T_2_m;
        rho_2_h = p_2_h / R_star / T_2_h;
        
        V_2A(end+1) = 3 * m / pi / b / D_m / (rho_2_t + rho_2_m + rho_2_h);
        end
        
        %   HOWELL CORRELATION %
        Re_How = 3e5;
        sigma_R_m = 1;
        
        %Chord calculation based on Howell value for Reynolds number
        