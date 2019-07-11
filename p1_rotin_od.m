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
        
        
    %%%% ROTOR INLET HIGH %%%%
    
    % Profiles
    delta_od_IGV(1) = delta_od_IGV_m;
    V_1T(1) =  V_1T_m;
    p_1(1) = p_1_m;
    
    alpha_1_geo_od = alpha_1_geo_high + IGV_rotation;
    
    
    s_IGV_high = pi * 2 * rhigh / N_IGV;
    
    o_IGV_high = s_IGV_high .* sind(alpha_1_geo_od);
    
    for i = 2:length(alpha_1_geo_od)
    
    T_T1(i) = T_T0;
      
                %p_1 CYCLE        
        p_1_iter = [p_1(i-1)*2 p_1(i-1)];
    while abs((p_1_iter(end)-p_1_iter(end-1))/p_1_iter(end-1)) > tol
        p_1_iter(end-1) = p_1_iter(end);
        
                       %Rho CYCLE
       rho_1_iter = [2*rho_1_m(end) rho_1_m(end)]; 
    
    while abs((rho_1_iter(end)-rho_1_iter(end-1))/rho_1_iter(end-1)) > tol
    
        rho_1_iter(end-1) = rho_1_iter(end); 
        
        T_1_iter = p_1_iter(end)/ R_star / rho_1_iter(end);
        V_1_iter = sqrt( 2*cp*( T_T1(i) - T_1_iter)  );  
        
            
       delta_od_IGV_iter = [2*delta_od_IGV_m 0*delta_od_IGV_m];
    while abs((delta_od_IGV_iter(end)-delta_od_IGV_iter(end-1))/delta_od_IGV_iter(end-1)) > tol 
        delta_od_IGV_iter(end-1) = delta_od_IGV_iter(end);
        
        alpha_1_iter = alpha_1_geo_od(i) - delta_od_IGV_iter(end);
        V_1A_iter = V_1_iter * cosd(alpha_1_iter);
        V_1T_iter = V_1_iter * sind(alpha_1_iter);
        W_1T_iter = V_1T_iter - U_t;
        W_1A_iter = V_1A_iter;
        W_1_iter = sqrt(W_1A_iter^2 + W_1T_iter^2);
        beta_1_iter = atand(W_1T_iter / W_1A_iter);
        
        delta_0_od_IGV_iter = asind(o_IGV_high(i)/s_IGV_high(i)*(1+(1-o_IGV_high(i)/s_IGV_high(i))*(alpha_1_geo_od(i)/90)^2)) - alpha_1_geo_od(i);
        
        Ma_1_iter_od = V_1_iter /sqrt(gamma*R_star*T_1_iter);
        
        if Ma_1_iter_od < 0.5
            delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter;
        else
            ics_delta_iter = 2 * Ma_1_iter_od - 1;
            delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter * (1 - 10*ics_delta_iter^3 + 15*ics_delta_iter^4 - 6*ics_delta_iter^5);
        end

    end
            
            V_1_iter = sqrt(V_1T_iter^2+V_1A_iter^2 );
            T_1_iter = T_T1(i) - V_1_iter^2 / 2/ cp;
            rho_1_iter(end+1) = p_1_iter(end) / R_star / T_1_iter;
  
    end
        
        
        Re_1_iter = V_1_iter * rho_1_iter(end) * c_IGV / mu;
 
        Y_1_p_tot_iter = y_AM_od(alpha_1_geo_od(i),alpha_0_t,alpha_1_iter,sigma_IGV_t,alpha_1_iter-alpha_1_t_geo,rho_1_iter(end),b,c_IGV_t,V_1_iter,mu);

    
    p_T1_iter = ( p_T0_t + Y_1_p_tot_iter * p_1_iter(end) ) / (Y_1_p_tot_iter+1);
    V_1_iter = sqrt( 2*gamma*R_star*T_1_iter / (gamma-1) * ( (p_T1_iter/p_1_iter(end))^((gamma-1)/gamma) - 1 ) );
            
            V_1T_iter = V_1_iter*sind(alpha_1_iter);
            
            T_1_iter = T_T1(i) - V_1_iter^2 / 2 / cp;
        
            p_1_iter(end+1) = p_1(i-1) + rho_1_iter(end) * V_1T_iter^2 / rhigh(i) *Dr;

    end
    
    delta_od_IGV(i) = delta_od_IGV_iter(end);
    V_1A(i) =  V_1A_iter;
    V_1T(i) =  V_1T_iter(end);
    V_1(i) =  V_1_iter(end);
    p_1(i) = p_1_iter(end);
    
    end

    V_1A_m(end+1) = (3*rho_0_m * V_0A - rho_1_iter(end)*V_1A_iter - rho_1_h(end)*V_1A_h)/rho_1_m(end);
    
    end
    

    


            story_rho_b   = b;
        story_rho_0   = rho_0;
        story_rho_1_m = rho_1_m;
        story_rho_1_t = rho_1_iter; 
        story_rho_0   = rho_0; 
        story_V_1A_m  = V_1A_m;
