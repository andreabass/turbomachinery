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