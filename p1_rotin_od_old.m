
%Assumption: we want to keep constant relative inlet angle at midspan

beta_1_m_od = beta_1_m;

V_1A_m = [2*V_1A(end) 174.68];
       
    while abs((V_1A_m(end)-V_1A_m(end-1))/V_1A_m(end-1)) > tol
       
    V_1A_m(end-1) = V_1A_m(end);  
    
    W_1A_m = V_1A_m(end);
    W_1T_m = W_1A_m * tand(beta_1_m_od);
    V_1T_m = W_1T_m + omega * D_m/2;
    alpha_1_m = atand(V_1T_m/V_1A_m(end));
    V_1_m = sqrt(V_1A_m(end)^2+V_1T_m^2);
    W_1_m = sqrt(W_1A_m^2+W_1T_m^2);
    
        
    %% MID
    
    T_T1_m = T_T0_m;
    T_1_m = T_T1_m - (V_1_m^2) / (2*cp); 
        
        %%%%%% [INITIALIZATION] %%%%%%%
        
        % IGV inlet density as initial value for rho_1_m
        
        rho_1_m = [2*rho_0(end) rho_0(end)]; 
        
    while abs((rho_1_m(end)-rho_1_m(end-1))/rho_1_m(end-1)) > tol
        
        rho_1_m(end-1) = rho_1_m(end); 
             
    Re_1_m = rho_1_m(end) * V_1_m * c_IGV / mu; 
    
    Y_1_p_tot_m = y_AM_od(alpha_0_m + (alpha_1_m - alpha_1_m_geo), alpha_0_m, alpha_1_m, sigma_IGV_m, alpha_1_m-alpha_1_m_geo, rho_1_m(end), b, c_IGV, V_1_m, mu);
    
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
%     alpha_1_m_geo_od = [2*alpha_1_m alpha_1_m];
%     while abs((alpha_1_m_geo_od(end)-alpha_1_m_geo_od(end-1))/alpha_1_m_geo_od(end-1)) > tol
%         alpha_1_m_geo_od(end-1) = alpha_1_m_geo_od(end);
%         
%         o_IGV_m = s_IGV_m * sind(alpha_1_m_geo_od(end));
%         
%         delta_0_od_IGV_m = asind(o_IGV_m/s_IGV_m*(1+(1-o_IGV_m/s_IGV_m)*(alpha_1_m_geo_od(end)/90)^2)) - alpha_1_m_geo_od(end);
%         
%         Ma_1_m_od = V_1_m /sqrt(gamma*R_star*T_1_m);
%         
%         if Ma_1_m_od < 0.5
%             delta_od_IGV_m = delta_0_od_IGV_m;
%         else
%             ics_delta = 2 * Ma_1_m_od - 1;
%             delta_od_IGV_m = delta_0_od_IGV_m * (1 - 10*ics_delta^3 + 15*ics_delta^4 - 6*ics_delta^5);
%         end
%         alpha_1_m_geo_od(end+1) = alpha_1_m - delta_od_IGV_m;
%     end
%         alpha_1_m_geo_od = alpha_1_m_geo_od(end);
         
%       IGV_rotation = alpha_1_m_geo_od - alpha_1_m_geo;

IGV_rotation = alpha_1_m - alpha_1_m_geo;
        
        
    %% ROTOR INLET high
    
    % Profiles
    delta_od_IGV_high = zeros(1,length(rhigh));
    %delta_od_IGV_high(1) = delta_od_IGV_m;
    V_1T_high(1) =  V_1T_m;
    p_1_high(1) = p_1_m;
    V_1A_high(1)    =  V_1A_m(end);
    V_1T_high(1)    =  V_1T_m(end);
    V_1_high(1)     =  V_1_m;
    alpha_1_high(1) = alpha_1_m;
    W_1A_high(1)    =  W_1A_m;
    W_1T_high(1)    =  W_1T_m;
    W_1_high(1)     =  W_1_m;
    beta_1_high(1)  = beta_1_m;
    p_1_high(1)     = p_1_m;
    T_1_high(1)     = T_1_m;
    rho_1_high(1)   = rho_1_m(end);
    
    alpha_1_geo_od_high = alpha_1_geo_high + IGV_rotation;
    
    s_IGV_high = pi * 2 * rhigh / N_IGV;
    
    sigma_IGV_high = c_IGV ./ s_IGV_high;
    
    o_IGV_high = s_IGV_high .* sind(alpha_1_geo_od_high);
    
    for i = 2:length(alpha_1_geo_od_high)
    
    T_T1_high(i) = T_T0;
    
%        delta_od_IGV_iter = [2*delta_od_IGV_high(i-1) delta_od_IGV_high(i-1)];
%     while abs((delta_od_IGV_iter(end)-delta_od_IGV_iter(end-1))/delta_od_IGV_iter(end-1)) > tol 
%         delta_od_IGV_iter(end-1) = delta_od_IGV_iter(end);
      
                %p_1 / T CYCLE        
        T_1_iter = [2*T_1_high(i-1) T_1_high(i-1)]; 
    while abs((T_1_iter(end)-T_1_iter(end-1))/T_1_iter(end-1)) > tol
        T_1_iter(end-1) = T_1_iter(end); 
            
        V_1_iter = sqrt( 2*cp*( T_T1_high(i) - T_1_iter(end) )  );
        
        alpha_1_iter = alpha_1_geo_od_high(i) - delta_od_IGV_high(i);
        % alpha_1_iter = alpha_1_geo_od_high(i) - delta_od_IGV_iter(end);
        V_1A_iter = V_1_iter * cosd(alpha_1_iter);
        V_1T_iter = V_1_iter * sind(alpha_1_iter);
        W_1T_iter = V_1T_iter - omega*rhigh(i);
        W_1A_iter = V_1A_iter;
        W_1_iter = sqrt(W_1A_iter^2 + W_1T_iter^2);
        beta_1_iter = atand(W_1T_iter / W_1A_iter);
    
       rho_1_iter = [rho_1_high(i-1)*2 rho_1_high(i-1)];
    while abs((rho_1_iter(end)-rho_1_iter(end-1))/rho_1_iter(end-1)) > tol
        rho_1_iter(end-1) = rho_1_iter(end);
        
        Re_1_iter = V_1_iter * rho_1_iter(end) * c_IGV / mu;
 
        Y_1_p_tot_iter = y_AM_od(alpha_0_m + IGV_rotation,alpha_0_t,alpha_1_iter,sigma_IGV_high(i),IGV_rotation,rho_1_iter(end),b,c_IGV,V_1_iter,mu);

        p_1_iter = p_1_high(i-1) + rho_1_iter(end) * V_1T_iter^2 / rhigh(i) * Drhigh;
        
        p_T1_iter = ( p_T0 + Y_1_p_tot_iter * p_1_iter ) / (Y_1_p_tot_iter+1);
        
        rho_1_iter(end+1) = p_1_iter / R_star / T_1_iter(end);     
    end  
    
    V_1_iter = sqrt( 2*gamma*R_star*T_1_iter(end) / (gamma-1) * ( (p_T1_iter/p_1_iter)^((gamma-1)/gamma) - 1 ) );
    
    T_1_iter(end+1) = T_T1_high(i) - V_1_iter^2/2/cp;
    
    end
    
%             delta_0_od_IGV_iter = asind(o_IGV_high(i)/s_IGV_high(i)*(1+(1-o_IGV_high(i)/s_IGV_high(i))*(alpha_1_geo_od_high(i)/90)^2)) - alpha_1_geo_od_high(i);
%         
%         Ma_1_iter_od = V_1_iter ./sqrt(gamma*R_star*T_1_iter(end))
%         
%         if Ma_1_iter_od < 0.5
%             delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter;
%         elseif Ma_1_iter_od >= 0.5 && Ma_1_iter_od < 1
%             ics_delta_iter = 2 * Ma_1_iter_od - 1;
%             delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter * (1 - 10*ics_delta_iter^3 + 15*ics_delta_iter^4 - 6*ics_delta_iter^5);
%         else
%             disp('Using 0° as supersonic deviation in high region')
%             delta_od_IGV_iter(end+1) = asind(o_IGV_high(i)/s_IGV_high(i)) - alpha_1_geo_od_high(i);
%         end
%         
%         
%     end
    
    % delta_od_IGV_high(i) = delta_od_IGV_iter(end);
    V_1A_high(i)    =  V_1A_iter;
    V_1T_high(i)    =  V_1T_iter(end);
    V_1_high(i)     =  V_1_iter(end);
    alpha_1_high(i) = alpha_1_iter(end);
    W_1A_high(i)    =  W_1A_iter;
    W_1T_high(i)    =  W_1T_iter(end);
    W_1_high(i)     =  W_1_iter(end);
    beta_1_high(i)  = beta_1_iter(end);
    p_1_high(i)     = p_1_iter(end);
    T_1_high(i)     = T_1_iter(end);
    rho_1_high(i)   = rho_1_iter(end);
    
    end
    
    %% ROTOR INLET LOW
    
    % Profiles
    delta_od_IGV_low = zeros(1,length(rlow));
    %delta_od_IGV_low(1) = delta_od_IGV_m;
    V_1T_low(1) =  V_1T_m;
    p_1_low(1) = p_1_m;
    V_1A_low(1)    =  V_1A_m(end);
    V_1T_low(1)    =  V_1T_m(end);
    V_1_low(1)     =  V_1_m;
    alpha_1_low(1) = alpha_1_m;
    W_1A_low(1)    =  W_1A_m;
    W_1T_low(1)    =  W_1T_m;
    W_1_low(1)     =  W_1_m;
    beta_1_low(1)  = beta_1_m;
    p_1_low(1)     = p_1_m;
    T_1_low(1)     = T_1_m;
    rho_1_low(1)   = rho_1_m(end);
    
    alpha_1_geo_od_low = alpha_1_geo_low + IGV_rotation;
    
    s_IGV_low = pi * 2 * rlow / N_IGV;
    
    sigma_IGV_low = c_IGV ./ s_IGV_low;
    
    o_IGV_low = s_IGV_low .* sind(alpha_1_geo_od_low);
    
    for i = 2:length(alpha_1_geo_od_low)
    
    T_T1_low(i) = T_T0;
    
%        delta_od_IGV_iter = [2*delta_od_IGV_low(i-1) delta_od_IGV_low(i-1)];
%     while abs((delta_od_IGV_iter(end)-delta_od_IGV_iter(end-1))/delta_od_IGV_iter(end-1)) > tol 
%         delta_od_IGV_iter(end-1) = delta_od_IGV_iter(end);
      
                %p_1 / T CYCLE        
        T_1_iter = [2*T_1_low(i-1) T_1_low(i-1)]; 
    while abs((T_1_iter(end)-T_1_iter(end-1))/T_1_iter(end-1)) > tol
        T_1_iter(end-1) = T_1_iter(end); 
            
        V_1_iter = sqrt( 2*cp*( T_T1_low(i) - T_1_iter(end))  );
        
        alpha_1_iter = alpha_1_geo_od_low(i) - delta_od_IGV_low(i);
        %alpha_1_iter = alpha_1_geo_od_low(i) - delta_od_IGV_iter(end);
        V_1A_iter = V_1_iter * cosd(alpha_1_iter);
        V_1T_iter = V_1_iter * sind(alpha_1_iter);
        W_1T_iter = V_1T_iter - omega*rlow(i);
        W_1A_iter = V_1A_iter;
        W_1_iter = sqrt(W_1A_iter^2 + W_1T_iter^2);
        beta_1_iter = atand(W_1T_iter / W_1A_iter);
    
       rho_1_iter = [rho_1_low(i-1)*2 rho_1_low(i-1)];
    while abs((rho_1_iter(end)-rho_1_iter(end-1))/rho_1_iter(end-1)) > tol
        rho_1_iter(end-1) = rho_1_iter(end);
        
        Re_1_iter = V_1_iter * rho_1_iter(end) * c_IGV / mu;
 
        Y_1_p_tot_iter = y_AM_od(alpha_0_m + IGV_rotation,alpha_0_t,alpha_1_iter,sigma_IGV_low(i),IGV_rotation,rho_1_iter(end),b,c_IGV,V_1_iter,mu);

        p_1_iter = p_1_low(i-1) + rho_1_iter(end) * V_1T_iter^2 / rlow(i) *Drlow;
        
        p_T1_iter = ( p_T0 + Y_1_p_tot_iter * p_1_iter ) / (Y_1_p_tot_iter+1);
    
        rho_1_iter(end+1) = p_1_iter / R_star / T_1_iter(end);     
    end  
    
    V_1_iter = sqrt( 2*gamma*R_star*T_1_iter(end) / (gamma-1) * ( (p_T1_iter/p_1_iter)^((gamma-1)/gamma) - 1 ) );
    
    T_1_iter(end+1) = T_T1_low(i) - V_1_iter^2/2/cp;
    
    end
    
%             delta_0_od_IGV_iter = asind(o_IGV_low(i)/s_IGV_low(i)*(1+(1-o_IGV_low(i)/s_IGV_low(i))*(alpha_1_geo_od_low(i)/90)^2)) - alpha_1_geo_od_low(i);
%         
%         Ma_1_iter_od = V_1_iter ./sqrt(gamma*R_star*T_1_iter(end))
%         
%         if Ma_1_iter_od < 0.5
%             delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter;
%         elseif Ma_1_iter_od >= 0.5 && Ma_1_iter_od < 1
%             ics_delta_iter = 2 * Ma_1_iter_od - 1;
%             delta_od_IGV_iter(end+1) = delta_0_od_IGV_iter * (1 - 10*ics_delta_iter^3 + 15*ics_delta_iter^4 - 6*ics_delta_iter^5);
%             disp('AA')
%         else
%             disp('Using 0° as supersonic deviation')
%             delta_od_IGV_iter(end+1) = asind(o_IGV_low(i)/s_IGV_low(i)) - alpha_1_geo_od_low(i);
%         end
%         
%         
%     end
    
    % delta_od_IGV_low(i) = delta_od_IGV_iter(end);
    
    V_1A_low(i)    =  V_1A_iter;
    V_1T_low(i)    =  V_1T_iter(end);
    V_1_low(i)     =  V_1_iter(end);
    alpha_1_low(i) = alpha_1_iter(end);
    W_1A_low(i)    =  W_1A_iter;
    W_1T_low(i)    =  W_1T_iter(end);
    W_1_low(i)     =  W_1_iter(end);
    beta_1_low(i)  = beta_1_iter(end);
    p_1_low(i)     = p_1_iter(end);
    T_1_low(i)     = T_1_iter(end);
    rho_1_low(i)   = rho_1_iter(end);
    
    end
    
    V_1A    =  [V_1A_low(end:-1:1) V_1A_high(2:end) ];
    V_1T    =  [V_1T_low(end:-1:1) V_1T_high(2:end) ];
    V_1     =  [V_1_low(end:-1:1) V_1_high(2:end) ];
    alpha_1 =  [alpha_1_low(end:-1:1) alpha_1_high(2:end) ];
    W_1A    =  [W_1A_low(end:-1:1) W_1A_high(2:end) ];
    W_1T    =  [W_1T_low(end:-1:1) W_1T_high(2:end) ];
    W_1     =  [W_1_low(end:-1:1) W_1_high(2:end) ];
    beta_1  =  [beta_1_low(end:-1:1) beta_1_high(2:end) ];
    p_1     =  [p_1_low(end:-1:1) p_1_high(2:end) ];
    T_1     =  [T_1_low(end:-1:1) T_1_high(2:end) ];
    rho_1   =  [rho_1_low(end:-1:1) rho_1_high(2:end) ];
    
%     V_1A_m(end+1) = ( m / pi / (b/2) - V_1A(4)*rho_1(4)*D_h/2 - V_1A(10)*rho_1(10)*D_t/2 ) / ( rho_1_m(end) * D_m );
    
      dA = 2 * pi * (r(1:end-1)+r(2:end))/2 .* Dr;
      rho_1_av = (rho_1(1:end-1)+rho_1(2:end))/2;
      V_1A_av = (V_1A(1:end-1)+V_1A(2:end))/2;
      dm = rho_1_av .* dA .* V_1A_av;
      mnew = sum(dm) - rho_1_m(end) .* ( pi * D_m .* Dr ) .* V_1A_m(end);
      dmnew_mid = m - mnew;
      V_1A_m(end+1) = dmnew_mid / ( rho_1_m(end) .* ( pi * D_m .* Dr ) );

    end
    


    


            story_rho_b   = b;
        story_rho_0   = rho_0;
        story_rho_1_m = rho_1_m;
        story_rho_1_t = rho_1_iter; 
        story_rho_0   = rho_0; 
        story_V_1A_m  = V_1A_m;
