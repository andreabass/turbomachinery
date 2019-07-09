%% ROTOR INLET - OFF DESIGN %%

%Assumptrion: we want to keep constant relative inlet angle at midspan

beta_1_m_od = beta_1_m;

V_1A_m = [2*V_1A(end) V_1A(end)];
       
    while abs((V_1A_m(end)-V_1A_m(end-1))/V_1A_m(end-1)) > tol
       
    V_1A_m(end-1) = V_1A_m(end);    
        
    %%%% MID %%%%
    
    W_1A_m = V_1A_m(end);
    W_1T_m = tand(beta_1_m_od) * W_1A_m;
    W_1_m = sqrt(W_1A_m^2 + W_1T_m^2);
    V_1T_m = U_m + W_1T_m;
    V_1_m = sqrt(V_1T_m(end)^2 + V_1A_m(end)^2);
    alpha_1_m = atand(V_1T_m(end) / V_1A_m(end));
    
    T_T1_m = T_T0_m;
    T_1_m = T_T1_m - (V_1_m^2) / (2*cp); 
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % IGV inlet density as initial value for rho_1_m
        
        rho_1_m = [2*rho_0(end) rho_0(end)]; 
        
    while abs((rho_1_m(end)-rho_1_m(end-1))/rho_1_m(end-1)) > tol
        
        rho_1_m(end-1) = rho_1_m(end); 
             
    Re_1_m = rho_1_m(end) * V_1_m * c_IGV / mu;
        
    Y_1_p_tot_m = y_AM_od(alpha_1_m_geo,alpha_0_m,alpha_1_m,sigma_IGV_m,alpha_1_m-alpha_1_m_geo,rho_1_m(end),b,c_IGV,V_1_m,mu);
    
    p_T1_m = [2*p_T0_m p_T0_m];
    while abs((p_T1_m(end)-p_T1_m(end-1))/p_T1_m(end-1))>tol
    p_T1_m(end-1) = p_T1_m(end);
    
    p_1_m = p_T1_m(end) / ( (1+(gamma-1)/2*(V_1_m^2/(gamma*R_star*T_1_m)))^(gamma/(gamma-1)) );
    
    p_T1_m(end+1) = ( p_T0_m + Y_1_p_tot_m * p_1_m ) / (Y_1_p_tot_m+1); 
    
    end
    p_T1_m = p_T1_m(end); 

    rho_1_m(end+1) = p_1_m / R_star / T_1_m;
    
    end
    
    %Cycle to determine delta, initializing the off design alpha angle
    alpha_1_m_geo_od = [2*(alpha_1_m + 3) alpha_1_m + 3];
    while abs((alpha_1_m_geo_od(end)-alpha_1_m_geo_od(end-1))/alpha_1_m_geo_od(end-1)) > tol
        alpha_1_m_geo_od(end-1) = alpha_1_m_geo_od(end);
        
        o_IGV_m = s_IGV_m * sind(alpha_1_m_geo_od(end));
        
        delta_0_od_IGV_m = asind(o_IGV_m/s_IGV_m*(1+(1-o_IGV_m/s_IGV_m)*(alpha_1_m_geo_od(end)/90)^2)) - alpha_1_m_geo_od(end);
        
        Ma_1_m_od = V_1_m /sqrt(gamma*R_star*T_1_m);
        
        if Ma_1_m_od < 0.5
            delta_od_IGV_m = delta_0_od_IGV_m;
        else
            ics_delta = 2 * Ma_1_m_od - 1;
            delta_od_IGV_m = delta_0_od_IGV_m * (1 - 10*ics_delta^3 + 15*ics_delta^4 - 6*ics_delta^5);
        end
        alpha_1_m_geo_od(end+1) = alpha_1_m - delta_od_IGV_m;
    end
        alpha_1_m_geo_od = alpha_1_m_geo_od(end);
        IGV_rotation = alpha_1_m_geo_od - alpha_1_m_geo;
        
        
    %%%% ROTOR INLET TIP %%%%
    alpha_1_t_geo_od = alpha_1_t_geo + IGV_rotation;
    
    T_T1_t = T_T0_t;
    
    delta_od_IGV_t = [2*delta_od_IGV_m delta_od_IGV_m];
    
    while abs((delta_od_IGV_t(end)-delta_od_IGV_t(end-1))/delta_od_IGV_t(end-1)) > tol
        delta_od_IGV_t(end-1) = delta_od_IGV_t(end);
        
        V_1T_t = [V_1T_m*2 V_1T_m];
        while abs((V_1T_t(end)-V_1T_t(end-1))/V_1T_t(end-1)) > tol
        V_1T_t(end-1) = V_1T_t(end);
        
        alpha_1_t = delta_od_IGV_t(end) + alpha_1_t_geo_od;
        V_1A_t = V_1T_t(end) / tand(alpha_1_t);
        V_1_t = sqrt(V_1A_t^2 + V_1T_t(end)^2);
        W_1T_t = V_1T_t(end) - U_t;
        W_1A_t = V_1A_t;
        W_1_t = sqrt(W_1A_t^2 + W_1T_t^2);
        beta_1_t = atand(W_1T_t / W_1A_t);
        
        %Rho initialization for losses calculation
        rho_1_t = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_t(end)-rho_1_t(end-1))/rho_1_t(end-1)) > tol
    
        rho_1_t(end-1) = rho_1_t(end); 
        
        Re_1_t = V_1_t * rho_1_t(end) * c_IGV / mu;
 
        Y_1_p_tot_t = y_AM_od(alpha_1_t_geo,alpha_0_t,alpha_1_t,sigma_IGV_t,alpha_1_t-alpha_1_t_geo,rho_1_t(end),b,c_IGV_t,V_1_t,mu);

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
        Ba = 1;
        V_1T_t(end+1) = sqrt((2/Ba*(p_1_t-p_1_m)/(D_t-D_m)*2 - (V_1T_m^2 * rho_1_m(end)/(D_m/2)))*D_t/2/rho_1_t(end));
        %V_1T_t(end+1) = V_1T_m * D_m/D_t;
    end
    
        o_IGV_t = s_IGV_t * sind(alpha_1_t_geo_od);
        
        delta_0_od_IGV_t = asind(o_IGV_t/s_IGV_t*(1+(1-o_IGV_t/s_IGV_t)*(alpha_1_t_geo_od/90)^2)) - alpha_1_t_geo_od;
        
        Ma_1_t_od = V_1_t /sqrt(gamma*R_star*T_1_t);
        
        if Ma_1_t_od < 0.5
            delta_od_IGV_t(end+1) = delta_0_od_IGV_t;
        else
            ics_delta_t = 2 * Ma_1_t_od - 1;
            delta_od_IGV_t(end+1) = delta_0_od_IGV_t * (1 - 10*ics_delta_t^3 + 15*ics_delta_t^4 - 6*ics_delta_t^5);
        end
    end

    %%%% ROTOR INLET HUB %%%%
    alpha_1_h_geo_od = alpha_1_h_geo + IGV_rotation;
    
    T_T1_h = T_T0_h;
    
    delta_od_IGV_h = [2*delta_od_IGV_m delta_od_IGV_m];
    
    while abs((delta_od_IGV_h(end)-delta_od_IGV_h(end-1))/delta_od_IGV_h(end-1)) > tol
        delta_od_IGV_h(end-1) = delta_od_IGV_h(end);
        
        V_1T_h = [V_1T_m*2 V_1T_m];
        while abs((V_1T_h(end)-V_1T_h(end-1))/V_1T_h(end-1)) > tol
        V_1T_h(end-1) = V_1T_h(end);
        
        alpha_1_h = delta_od_IGV_h(end) + alpha_1_h_geo_od;
        V_1A_h = V_1T_h(end) / tand(alpha_1_h);
        V_1_h = sqrt(V_1A_h^2 + V_1T_h(end)^2);
        W_1T_h = V_1T_h(end) - U_h;
        W_1A_h = V_1A_h;
        W_1_h = sqrt(W_1A_h^2 + W_1T_h^2);
        beta_1_h = atand(W_1T_h / W_1A_h);
        
        %Rho initialization for losses calculation
        rho_1_h = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_h(end)-rho_1_h(end-1))/rho_1_h(end-1)) > tol
    
        rho_1_h(end-1) = rho_1_h(end); 
        
        Re_1_h = V_1_h * rho_1_h(end) * c_IGV_h / mu;
 
        Y_1_p_tot_h = y_AM_od(alpha_1_h_geo,alpha_0_h,alpha_1_h,sigma_IGV_h,alpha_1_h-alpha_1_h_geo,rho_1_h(end),b,c_IGV_h,V_1_h,mu);

        T_1_h = T_T1_h - (V_1_h^2) / (2*cp);
    
    p_T1_h = [2*p_T0_h p_T0_h];
    while abs((p_T1_h(end)-p_T1_h(end-1))/p_T1_h(end-1))>tol
    p_T1_h(end-1) = p_T1_h(end);
    
    p_1_h = p_T1_h(end) / ( (1+(gamma-1)/2*(V_1_h^2/(gamma*R_star*T_1_h)))^(gamma/(gamma-1)) );
    
    p_T1_h(end+1) = ( p_T0_h + Y_1_p_tot_h * p_1_h ) / (Y_1_p_tot_h+1);
    
    end
    p_T1_h = p_T1_h(end);
    
    rho_1_h(end+1) = p_1_h / R_star / T_1_h;
    
    end
        Ba = 1;
        V_1T_h(end+1) = sqrt((2/Ba* (p_1_h-p_1_m)/(D_h-D_m) *2 - (V_1T_m^2 * rho_1_m(end)/(D_m/2)))*D_h/2/rho_1_h(end));
        %V_1T_h(end+1) = V_1T_m * D_m/D_h;
    end
    
        o_IGV_h = s_IGV_h * sind(alpha_1_h_geo_od);
        
        delta_0_od_IGV_h = asind(o_IGV_h/s_IGV_h*(1+(1-o_IGV_h/s_IGV_h)*(alpha_1_h_geo_od/90)^2)) - alpha_1_h_geo_od;
        
        Ma_1_h_od = V_1_h /sqrt(gamma*R_star*T_1_h);
        
        if Ma_1_h_od < 0.5
            delta_od_IGV_h(end+1) = delta_0_od_IGV_h;
        else
            ics_delta_h = 2 * Ma_1_t_od - 1;
            delta_od_IGV_h(end+1) = delta_0_od_IGV_h * (1 - 10*ics_delta_h^3 + 15*ics_delta_h^4 - 6*ics_delta_h^5);
        end
    end
    
    V_1A_m(end+1) = (3*rho_0_m * V_0A - rho_1_t(end)*V_1A_t - rho_1_h(end)*V_1A_h)/rho_1_m(end);
    
    end

            story_rho_b   = b;
        story_rho_0   = rho_0;
        story_rho_1_m = rho_1_m;
        story_rho_1_t = rho_1_t; 
        story_rho_0   = rho_0; 
        story_V_1A_m  = V_1A_m;
