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
         
        V_1A = [V_0A*2 V_0A*0.5];
       
    while abs((V_1A(end)-V_1A(end-1))/V_1A(end-1)) > tol    
            V_1A(end-1) = V_1A(end); 
        
         T_1_t = [300*2 300];
       
    while abs((T_1_t(end)-T_1_t(end-1))/T_1_t(end-1)) > tol    
            T_1_t(end-1) = T_1_t(end); 
            
       W_1_t =  M_R1_t * sqrt(gamma * R_star * T_1_t(end));
       W_1A_t = V_1A(end);
       W_1T_t = - sqrt(W_1_t^2-W_1A_t^2);
       V_1T_t = W_1T_t + omega*D_t/2;
       V_1_t = sqrt(V_1A(end)^2 + V_1T_t^2);
       T_T1_t = T_T0_t;
       T_1_t(end+1) = T_T1_t - V_1_t^2/2/cp;
       
    end
    T_1_t = T_1_t(end);
    
       alpha_1_t = tand(V_1T_t / V_1A(end));
       beta_1_t = tand(W_1T_t/W_1A_t);
      
       V_1T_h = V_1T_t * D_t / D_h;
       V_1_h = sqrt(V_1T_h^2 + V_1A(end)^2);
       W_1T_h = V_1T_h - omega*D_h/2;
       W_1A_h = V_1A(end);
       W_1_h = sqrt(W_1T_h^2 + W_1A_h^2);
       alpha_1_h = tand(V_1T_h/V_1A(end));
       beta_1_h = tand(W_1T_h/W_1A_h); 
       
       V_1T_m = V_1T_t * D_t / D_m;
       V_1_m = sqrt(V_1T_m^2 + V_1A(end)^2);
       W_1T_m = V_1T_m - omega*D_m/2;
       W_1A_m = V_1A(end);
       W_1_m = sqrt(W_1T_m^2 + W_1A_m^2);
       alpha_1_m = tand(V_1T_m/V_1A(end)); 
       beta_1_m = tand(W_1T_m/W_1A_m); 

 
             
    Y_p_1_in_m = y_AM_inc_min(alpha_1_m);
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % IGV inlet density as initial value for rho_1_m
        
        rho_1_m = [2*rho_0(end) rho_0(end)]; 
        
    while abs((rho_1_m(end)-rho_1_m(end-1))/rho_1_m(end-1)) > tol
        
        rho_1_m(end-1) = rho_1_m(end); 
             
    Re_1_m = rho_1_m(end) * V_1_m * c_IGV / mu;
    
    Y_p_1_Re_m = y_AM_Re(Y_p_1_in_m,Re_1_m); 
    
    s_over_c_min_m = s_c_min_AM(alpha_1_m);

    Y_1_sec_m = y_AM_sec(alpha_1_m,alpha_0_m,c_IGV,b,s_over_c_min_m);
        
    Y_1_p_tot_m = Y_p_1_Re_m + Y_1_sec_m;
    
    T_T1_m = T_T0_m;
    
    T_1_m = T_T1_m - (V_1_m^2) / (2*cp);
    
    p_T1_m = [2*p_T0_m p_T0_m];
    while abs((p_T1_m(end)-p_T1_m(end-1))/p_T1_m(end-1))>tol
    p_T1_m(end-1) = p_T1_m(end);
    
    p_1_m = p_T1_m(end) / ( (1+(gamma-1)/2*(V_1_m^2/(gamma*R_star*T_1_m)))^(gamma/(gamma-1)) );
    
    p_T1_m(end+1) = ( p_T0_m + Y_1_p_tot_m * p_1_m ) / (Y_1_p_tot_m+1); 
    
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

    Y_p_1_in_t = y_AM_inc(alpha_1_t,s_over_c_t);

        %%%%%% [INITIALIZATION] %%%%%%%
        
        % Current density @ midspan as initial value for rho_1_t
        
        rho_1_t = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_t(end)-rho_1_t(end-1))/rho_1_t(end-1)) > tol
    
        rho_1_t(end-1) = rho_1_t(end); 
        
    Re_1_t = V_1_t * rho_1_t(end) * c_IGV / mu;
    
    Y_p_1_Re_t = y_AM_Re(Y_p_1_in_t,Re_1_t);
    
    Y_1_sec_t = y_AM_sec(alpha_1_t,alpha_0_t,c_IGV,b,s_over_c_t);
    
    Y_1_p_tot_t = Y_p_1_Re_t + Y_1_sec_t;
    
    T_T1_t = T_T0_t;

    T_1_t = T_T1_t - (V_1_t^2) / (2*cp);
    
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

    Y_p_1_in_h = y_AM_inc(alpha_1_h,s_over_c_h);
    
    T_T1_h = T_T0_h;

    T_1_h = T_T1_h - (V_1_h^2) / (2*cp);
    
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % Current density @ midspan as initial value for rho_1_t
        
        rho_1_h = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_h(end)-rho_1_h(end-1))/rho_1_h(end-1)) > tol
    
        rho_1_h(end-1) = rho_1_h(end); 
    
    Re_1_h = rho_1_h(end) * c_IGV * V_1_h / mu;
        
    Y_p_1_Re_h = y_AM_Re(Y_p_1_in_h,Re_1_h);
    
    Y_1_sec_h = y_AM_sec(alpha_1_h,alpha_0_h,c_IGV,b,s_over_c_h);
    
    Y_1_p_tot_h = Y_p_1_Re_h + Y_1_sec_h;
    
    p_T1_h = [2*p_T0_h p_T0_h];
    while abs((p_T1_h(end)-p_T1_h(end-1))/p_T1_h(end-1))>tol
    p_T1_h(end-1) = p_T1_h(end);
    
    p_1_h = p_T1_h(end) / ( (1+(gamma-1)/2*(V_1_h^2/(gamma*R_star*T_1_h)))^(gamma/(gamma-1)) );
    
    p_T1_h(end+1) = ( p_T0_h + Y_1_p_tot_h * p_1_h ) / (Y_1_p_tot_h+1);  
    
    end
    p_T1_h = p_T1_h(end);
    
    rho_1_h(end+1) = p_1_h / R_star / T_1_h;
    
    end

    V_1A(end+1) = 3 * rho_0(end) * V_0A / ( rho_1_h(end) + rho_1_m(end) + rho_1_t(end) );
    
    end
        
        
        V_1A = V_1A(end);
        V_1A_m = V_1A;
        V_1A_t = V_1A;
        V_1A_h = V_1A;
        rho_0   = rho_0(end);
        rho_1_m = rho_1_m(end);
        rho_1_t = rho_1_t(end);