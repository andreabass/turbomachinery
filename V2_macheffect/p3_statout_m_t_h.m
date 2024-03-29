         %% STATOR  OUTLET %%
          %%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%
             %%%%%%%%%%%
              %%%%%%%%%
               %%%%%%%
                %%%%%
                 %%%
                  %
              
        % Velocity triangle @ stator outlet
        
        V_3T_m = V2Tm_V1Tm * V_2T_m;
        
        V_3T_t = V_3T_m * D_m / D_t;
        V_3T_h = V_3T_m * D_m / D_h(end);
        
        %%% INITIALIZATION %%%

        % We initialize the losses in the rotor by assuming each static
        % efficiency (mid, tip. hub) equal to the total-to-total efficiency
        % of the stage
        
        eta_S_t = [eta_TT(end) eta_TT(end) + 2*tol];
        eta_S_m = [eta_TT(end) eta_TT(end) + 2*tol];
        eta_S_h = [eta_TT(end) eta_TT(end) + 2*tol];    
        
        % Constant axial velocity as initial value
        
        V_3A = V_2A;
        V_3A = [V_3A V_3A + 2*tol];
        
        for i=1:KK
        while abs(eta_S_t(end) - eta_S_t(end-1))>tol || abs(eta_S_m(end) - eta_S_m(end-1))>tol || abs(eta_S_h(end) - eta_S_h(end-1))>tol
            eta_S_t(end-1) = eta_S_t(end);
            eta_S_m(end-1) = eta_S_m(end);
            eta_S_h(end-1) = eta_S_h(end);
        
            p3_statvelout
        
        V_3A_t = V_3A;
        V_3A_h = V_3A;
        
        M_2_m = V_2_m / sqrt(gamma*R_star*T_2_m);
        M_2_h = V_2_h / sqrt(gamma*R_star*T_2_h);
        M_2_t = V_2_t / sqrt(gamma*R_star*T_2_t);
        
        %   HOWELL CORRELATION %
        %These assumptions are used to obtain the two correction
        %coefficients of Howell correlation (Psi, Phi) equal to one. In
        %principle is necessary to check also that the deflection Dbeta is
        %sufficiently close to the optimal value associated to beta2
            
        Re_How = 3e5;
  
           if i==1
             sigma_S_m = sigma_S_m_design;              
           end
        
        %Chord calculation based on Howell value for Reynolds number
        
            c_How_m_S = mu*Re_How/rho_3_m/V_3_m;
            c_How_t_S = mu*Re_How/rho_3_t/V_3_t;
            c_How_h_S = mu*Re_How/rho_3_h/V_3_h;
         c_S_m = max( max(max(c_How_m_S,c_How_t_S), c_How_h_S), c_S_design );
         c_S_t = c_S_m;
         c_S_h = c_S_m;

         %Rotor geometry
         
         s_S_m = c_S_m / sigma_S_m;
         N_S = ceil(pi * D_m / s_S_m);
            s_S_m = pi * D_m / N_S;
            sigma_S_m = c_S_m / s_S_m;
          
        s_S_t = pi * D_t / N_S;
        sigma_S_t = c_S_t / s_S_t;
        s_S_h = pi * D_h(end) / N_S;
        sigma_S_h = c_S_h / s_S_h;
        
        geo_stator
        
        % Stall angles (MID)
        
    attack_S_c_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_c_m(end) - attack_S_c_m(end-1))/attack_S_c_m(end-1) ) > tol
        attack_S_c_m(end-1) = attack_S_c_m(end);
        
        alpha_c_m = -attack_S_c_m(end) + gamma_statm;
        
        if abs(alpha_c_m)<20
            alpha_c_m = 20;
        end
        
        attack_S_c_m(end+1) = attack_S_m_design - 9 + ( 1 - (30/abs(alpha_c_m))^0.48 )*teta_mid_stat/4.176;
        
    end
    attack_S_c_m = attack_S_c_m(end);
    
    R_S_cm = attack_S_m_design - attack_S_c_m;
    
    i_S_c_m = i_opt_statm - R_S_cm / (1+0.5*M_2_m^3);
    
   attack_S_s_m = [attack_S_m_design*2 attack_S_m_design];
    while abs( (attack_S_s_m(end) - attack_S_s_m(end-1))/attack_S_s_m(end-1) ) > tol
        attack_S_s_m(end-1) = attack_S_s_m(end);
        
        alpha_s_m = -attack_S_s_m(end) + gamma_statm;
        attack_S_s_m(end+1) = attack_S_m_design + 10.3 + ( 2.92 - (abs(alpha_s_m)/15.6) )*teta_mid_stat/8.2;
        
    end
    attack_S_s_m = attack_S_s_m(end);
    
    R_S_sm = - attack_S_m_design + attack_S_s_m;
    
    i_S_s_m = i_opt_statm + R_S_sm / (1+0.5*(Ksh_i*M_2_m)^3);
    
    i_S_mm  = i_S_c_m + (i_S_s_m - i_S_c_m)*R_S_cm/(R_S_cm+R_S_sm);
    
        % Stall angles (HUB)
        
    attack_S_c_h = [attack_S_h_design*2 attack_S_h_design];
    while abs( (attack_S_c_h(end) - attack_S_c_h(end-1))/attack_S_c_h(end-1) ) > tol
        attack_S_c_h(end-1) = attack_S_c_h(end);

        beta_S_c_h = -attack_S_c_h(end) + gamma_stath;
        
                if abs(beta_S_c_h)<20
            beta_S_c_h = 20;
                end
        
        attack_S_c_h(end+1) = attack_S_h_design - 9 + ( 1 - (30/abs(beta_S_c_h))^0.48 )*teta_hub_stat/4.176;
        
    end
    attack_S_c_h = attack_S_c_h(end);
    
    R_S_ch = attack_S_h_design - attack_S_c_h;
    
    i_S_c_h = i_opt_stath - R_S_ch / (1+0.5*M_2_h^3);
    
   attack_S_s_h = [attack_S_h_design*2 attack_S_h_design];
    while abs( (attack_S_s_h(end) - attack_S_s_h(end-1))/attack_S_s_h(end-1) ) > tol
        attack_S_s_h(end-1) = attack_S_s_h(end);
        
        alpha_s_h = -attack_S_s_h(end) + gamma_stath;
        attack_S_s_h(end+1) = attack_S_h_design + 10.3 + ( 2.92 - (abs(alpha_s_h)/15.6) )*teta_hub_stat/8.2;
        
    end
    attack_S_s_h = attack_S_s_h(end);
    
    R_S_sh = - attack_S_h_design + attack_S_s_h;
    
    i_S_s_h = i_opt_stath + R_S_sh / (1+0.5*(Ksh_i*M_2_h)^3);
    
    i_S_mh  = i_S_c_h + (i_S_s_h - i_S_c_h)*R_S_ch/(R_S_ch+R_S_sh);
    
            % Stall angles (TIP)
        
    attack_S_c_t = [attack_S_h_design*2 attack_S_h_design];
    while abs( (attack_S_c_t(end) - attack_S_c_t(end-1))/attack_S_c_t(end-1) ) > tol
        attack_S_c_t(end-1) = attack_S_c_t(end);

        alpha_c_t = -attack_S_c_t(end) + gamma_statt;
        
                if abs(alpha_c_t)<20
            alpha_c_t = 20;
                end
        
        attack_S_c_t(end+1) = attack_S_t_design - 9 + ( 1 - (30/abs(alpha_c_t))^0.48 )*teta_tip_stat/4.176;
        
    end
    attack_S_c_t = attack_S_c_t(end);
    
    R_S_ct = attack_S_t_design - attack_S_c_t;
    
    i_S_c_t = i_opt_statt - R_S_ct / (1+0.5*M_2_t^3);
    
   attack_S_s_t = [attack_S_t_design*2 attack_S_t_design];
    while abs( (attack_S_s_t(end) - attack_S_s_t(end-1))/attack_S_s_t(end-1) ) > tol
        attack_S_s_t(end-1) = attack_S_s_t(end);
        
        alpha_s_t = -attack_S_s_t(end) + gamma_statt;
        attack_S_s_t(end+1) = attack_S_t_design + 10.3 + ( 2.92 - (abs(alpha_s_t)/15.6) )*teta_tip_stat/8.2;
        
    end
    attack_S_s_t = attack_S_s_t(end);
    
    R_S_st = - attack_S_t_design + attack_S_s_t;
    
    i_S_s_t = i_opt_statt + R_S_st / (1+0.5*(Ksh_i*M_2_t)^3);
    
    i_S_mt  = i_S_c_t + (i_S_s_t - i_S_c_t)*R_S_ct/(R_S_ct+R_S_st);
    
        %Losses correlations
        
        K_S_m = abs(tand(alpha_2_m) - V_3A(end)/V_2A_m * tand(alpha_3_m));
        K_S_t = abs(tand(alpha_2_t) - V_3A(end)/V_2A_t * tand(alpha_3_t));
        K_S_h = abs(tand(alpha_2_h) - V_3A(end)/V_2A_h * tand(alpha_3_h));
        
        Vmax_V2_m = 1.12 + 0.61 * cosd(alpha_2_m)^2/sigma_S_m * K_S_m;
        Vmax_V2_t = 1.12 + 0.61 * cosd(alpha_2_t)^2/sigma_S_t * K_S_t;
        Vmax_V2_h = 1.12 + 0.61 * cosd(alpha_2_h)^2/sigma_S_h * K_S_h;
        
        DmS = Vmax_V2_m * V_2_m / V_3_m;
        DtS = Vmax_V2_t * V_2_t / V_3_t;
        DhS = Vmax_V2_h * V_2_h / V_3_h;
        
        Y_3_p_tot_m = 0.004 * ( 1 + 3.1*(DmS-1)^2 + 0.4*(DmS-1)^8 ) * 2 * sigma_S_m / cosd(alpha_3_m) * (V_3_m/V_2_m)^2;
        Y_3_p_tot_t = 0.004 * ( 1 + 3.1*(DtS-1)^2 + 0.4*(DtS-1)^8 ) * 2 * sigma_S_t / cosd(alpha_3_t) * (V_3_t/V_2_t)^2;
        Y_3_p_tot_h = 0.004 * ( 1 + 3.1*(DhS-1)^2 + 0.4*(DhS-1)^8 ) * 2 * sigma_S_h / cosd(alpha_3_h) * (V_3_h/V_2_h)^2;
        
        Y_3_p_tot_m = Y_3_p_tot_m * ( 1 + (i_S_mm - i_opt_statm)^2 / R_S_sm^2 ); 
        Y_3_p_tot_h = Y_3_p_tot_h * ( 1 + (i_S_mh - i_opt_stath)^2 / R_S_sh^2 );
        Y_3_p_tot_t = Y_3_p_tot_t * ( 1 + (i_S_mt - i_opt_statt)^2 / R_S_st^2 );
        
        % Critical Mach calculation
        
        T_cr_2_m = T_T2_m * 2 / (gamma+1);
        T_cr_2_h = T_T2_h * 2 / (gamma+1);
        T_cr_2_t = T_T2_t * 2 / (gamma+1);
        
        V_cr_2_m   = 1 * sqrt(gamma*R_star*T_cr_2_m);
        V_cr_2_h   = 1 * sqrt(gamma*R_star*T_cr_2_h);
        V_cr_2_t   = 1 * sqrt(gamma*R_star*T_cr_2_t);
        
        M_cr_2_m   = M_2_m * V_cr_2_m / ( V_2_m * Vmax_V2_m );
        M_cr_2_h   = M_2_h * V_cr_2_h / ( V_2_h * Vmax_V2_h );
        M_cr_2_t   = M_2_t * V_cr_2_t / ( V_2_t * Vmax_V2_t );
        
        if M_2_m > M_cr_2_m
        Y_3_p_tot_m = Y_3_p_tot_m + Ksh_i * ((M_2_m/M_cr_2_m-1)*V_cr_2_m/V_2_m)^2;
        end
        
        if M_2_h > M_cr_2_h
        Y_3_p_tot_h = Y_3_p_tot_h + Ksh_i * ((M_2_h/M_cr_2_h-1)*V_cr_2_h/V_2_h)^2;
        end
        
        if M_2_t > M_cr_2_t
        Y_3_p_tot_t = Y_3_p_tot_t + Ksh_i * ((M_2_t/M_cr_2_t-1)*V_cr_2_t/V_2_t)^2;
        end
        
        %TDN variables
        
        p_T3_m = p_T2_m - Y_3_p_tot_m * (p_T2_m - p_2_m);
        p_T3_t = p_T2_t - Y_3_p_tot_t * (p_T2_t - p_2_t);
        p_T3_h = p_T2_h - Y_3_p_tot_h * (p_T2_h - p_2_h);  
        
        p_3_m = p_T3_m / ( (1+(gamma-1)/2*(rho_3_m*V_3_m^2/(gamma*p_3_m)))^(gamma/(gamma-1)) );
        p_3_t = p_T3_t / ( (1+(gamma-1)/2*(rho_3_t*V_3_t^2/(gamma*p_3_t)))^(gamma/(gamma-1)) );
        p_3_h = p_T3_h / ( (1+(gamma-1)/2*(rho_3_h*V_3_h^2/(gamma*p_3_h)))^(gamma/(gamma-1)) );

        T_3_m_is  = T_2_m * ( p_3_m / p_2_m )^((gamma-1)/gamma);
        T_3_t_is  = T_2_t * ( p_3_t / p_2_t )^((gamma-1)/gamma);
        T_3_h_is  = T_2_h * ( p_3_h / p_2_h )^((gamma-1)/gamma);
        
        T_T3_m_is  = T_T2_m * ( p_T3_m / p_T2_m )^((gamma-1)/gamma);
        T_T3_t_is  = T_T2_t * ( p_T3_t / p_T2_t )^((gamma-1)/gamma);
        T_T3_h_is  = T_T2_h * ( p_T3_h / p_T2_h )^((gamma-1)/gamma);
        
        %Update SS efficiency @ rotor outlet
        
        eta_S_m(end+1) = ( T_3_m_is - T_2_m ) / ( T_3_m - T_2_m );
        eta_S_t(end+1) = ( T_3_t_is - T_2_t ) / ( T_3_t - T_2_t );
        eta_S_h(end+1) = ( T_3_h_is - T_2_h ) / ( T_3_h - T_2_h );
        
        %Update velocity triangle @ rotor outlet
        
            rho_3_m = p_3_m / R_star / T_3_m; 
            rho_3_t = p_3_t / R_star / T_3_t; 
            rho_3_h = p_3_h / R_star / T_3_h; 
            
            V_3A(end+1) = 3 * m / pi / b / D_m / (rho_3_t + rho_3_m + rho_3_h);        

        end
        
        %Now we optimize the solidity with the Howell correlation and we
        %evaluate again the efficiency
        if i==1 && HOW_OPT == 1
           Db_Psi = ppval(Dbeta_Psi_curve, abs(alpha_3_m));
           Psi_opt = abs(alpha_3_m-alpha_2_m)/Db_Psi;          
           if Psi_opt<1.329 && Psi_opt>0.749
           x=0.4:0.001:1.6;
           s_over_c_S = mean(x( find(ppval(Psi_curve, x)>Psi_opt-0.001 & ppval(Psi_curve, x)<Psi_opt+0.001)));
           Howell_S='Stator Optimized';
           elseif Psi_opt>1.329
           s_over_c_S = 0.4;
           Howell_S='Stator Overloaded';
           else
           s_over_c_S = 1.6;        
           Howell_S='Stator Underloaded';
           end
           sigma_S_m = 1/s_over_c_S;
           eta_S_t = [eta_S_t(end) eta_S_t(end) + 2*tol];
           eta_S_m = [eta_S_m(end) eta_S_m(end) + 2*tol];
           eta_S_h = [eta_S_h(end) eta_S_h(end) + 2*tol];   
        end
        
        end
        
        
        story_eta_S_t = eta_S_t;
        story_eta_S_m = eta_S_m; 
        story_eta_S_h = eta_S_h; 
        story_V_3A    = V_3A;

        eta_S_m   = eta_S_m(end);
        eta_S_t   = eta_S_t(end);
        eta_S_h   = eta_S_h(end);
        V_3      = V_3A(end);
        
        
        
        
           

