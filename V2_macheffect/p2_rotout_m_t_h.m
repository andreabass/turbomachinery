         %% ROTOR  OUTLET %%
          %%%%%%%%%%%%%%%%%
           %%%%%%%%%%%%%%%
            %%%%%%%%%%%%%
             %%%%%%%%%%%
              %%%%%%%%%
               %%%%%%%
                %%%%%
                 %%%
                  %
                  
        l_Eu = deltaHis_TT / eta_TT(end);
                     
        %%% INITIALIZATION %%%

        % We initialize the losses in the rotor by assuming each static
        % efficiency (mid, tip. hub) equal to the total-to-total efficiency
        % of the stage
        
        eta_R_t = [2*eta_TT(end) eta_TT(end)];
        eta_R_m = [2*eta_TT(end) eta_TT(end)];
        eta_R_h = [2*eta_TT(end) eta_TT(end)];    
        
        % Constant axial velocity as initial value
        
        V_2A = V_1A_m;
        V_2A = [V_2A V_2A + 2*tol];
        Howell_correlation
        opt=1;
        
        % Check if the Howell optimization procedure was selected
        if HOW_OPT == 1
            KK = 2;
        else
            KK = 1;
        end
        
        for i=1:KK
        while abs(eta_R_t(end) - eta_R_t(end-1))>tol || abs(eta_R_m(end) - eta_R_m(end-1))>tol || abs(eta_R_h(end) - eta_R_h(end-1))>tol
            eta_R_t(end-1) = eta_R_t(end);
            eta_R_m(end-1) = eta_R_m(end);
            eta_R_h(end-1) = eta_R_h(end);
        
            p2_rotvelout
            
        M_R1_m = W_1_m / sqrt(gamma*R_star*T_1_m);
        M_R1_h = W_1_h / sqrt(gamma*R_star*T_1_h);
        M_R1_t = W_1_t / sqrt(gamma*R_star*T_1_t);
        
        %   HOWELL CORRELATION %
        %These assumptions are used to obtain the two correction
        %coefficients of Howell correlation (Psi, Phi) equal to one. In
        %principle is necessary to check also that the deflection Dbeta is
        %sufficiently close to the optimal value associated to beta2
            
        Re_How = 3e5;
       
        if i==1
            sigma_R_m = sigma_R_m_design; 
        end
        
        
        %Chord calculation based on Howell value for Reynolds number
        
            c_How_m = mu*Re_How/rho_2_m/W_2_m;
            c_How_t = mu*Re_How/rho_2_t/W_2_t;
            c_How_h = mu*Re_How/rho_2_h/W_2_h;
         c_R_m = max( max(max(c_How_m,c_How_t), c_How_h), c_R_design );
         c_R_t = c_R_m;
         c_R_h = c_R_m;

         %Rotor geometry
         
         s_R_m = c_R_m / sigma_R_m;
         N_R = ceil(pi * D_m / s_R_m);
            s_R_m = pi * D_m / N_R;
            sigma_R_m = c_R_m / s_R_m;
          
        s_R_t = pi * D_t / N_R;
        sigma_R_t = c_R_t / s_R_t;
        s_R_h = pi * D_h(end) / N_R;
        sigma_R_h = c_R_h / s_R_h;
        
        geo_rotor
        
        % Stall angles (MID)
        
    attack_c_m = [attack_m_design*2 attack_m_design];
    while abs( (attack_c_m(end) - attack_c_m(end-1))/attack_c_m(end-1) ) > tol
        attack_c_m(end-1) = attack_c_m(end);
        
        beta_c_m = -attack_c_m(end) + gamma_rotm;
        
        if abs(beta_c_m)<20
            beta_c_m = 20;
        end
        
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
    
    i_s_m = i_opt_rotm + R_sm / (1+0.5*(Ksh_i*M_R1_m)^3);
    
    i_mm  = i_c_m + (i_s_m - i_c_m)*R_cm/(R_cm+R_sm);
    
        % Stall angles (HUB)
        
    attack_c_h = [attack_h_design*2 attack_h_design];
    while abs( (attack_c_h(end) - attack_c_h(end-1))/attack_c_h(end-1) ) > tol
        attack_c_h(end-1) = attack_c_h(end);

        beta_c_h = -attack_c_h(end) + gamma_roth;
        
                if abs(beta_c_h)<20
            beta_c_h = 20;
                end
        
        attack_c_h(end+1) = attack_h_design - 9 + ( 1 - (30/abs(beta_c_h))^0.48 )*teta_hub/4.176;
        
    end
    attack_c_h = attack_c_h(end);
    
    R_ch = attack_h_design - attack_c_h;
    
    i_c_h = i_opt_roth - R_ch / (1+0.5*M_R1_h^3);
    
   attack_s_h = [attack_h_design*2 attack_h_design];
    while abs( (attack_s_h(end) - attack_s_h(end-1))/attack_s_h(end-1) ) > tol
        attack_s_h(end-1) = attack_s_h(end);
        
        beta_s_h = -attack_s_h(end) + gamma_roth;
        attack_s_h(end+1) = attack_h_design + 10.3 + ( 2.92 - (abs(beta_s_h)/15.6) )*teta_hub/8.2;
        
    end
    attack_s_h = attack_s_h(end);
    
    R_sh = - attack_h_design + attack_s_h;
    
    i_s_h = i_opt_roth + R_sh / (1+0.5*(Ksh_i*M_R1_h)^3);
    
    i_mh  = i_c_h + (i_s_h - i_c_h)*R_ch/(R_ch+R_sh);
    
            % Stall angles (TIP)
        
    attack_c_t = [attack_h_design*2 attack_h_design];
    while abs( (attack_c_t(end) - attack_c_t(end-1))/attack_c_t(end-1) ) > tol
        attack_c_t(end-1) = attack_c_t(end);

        beta_c_t = -attack_c_t(end) + gamma_rott;
        
                if abs(beta_c_t)<20
            beta_c_t = 20;
                end
        
        attack_c_t(end+1) = attack_t_design - 9 + ( 1 - (30/abs(beta_c_t))^0.48 )*teta_tip/4.176;
        
    end
    attack_c_t = attack_c_t(end);
    
    R_ct = attack_t_design - attack_c_t;
    
    i_c_t = i_opt_rott - R_ct / (1+0.5*M_R1_t^3);
    
   attack_s_t = [attack_t_design*2 attack_t_design];
    while abs( (attack_s_t(end) - attack_s_t(end-1))/attack_s_t(end-1) ) > tol
        attack_s_t(end-1) = attack_s_t(end);
        
        beta_s_t = -attack_s_t(end) + gamma_rott;
        attack_s_t(end+1) = attack_t_design + 10.3 + ( 2.92 - (abs(beta_s_t)/15.6) )*teta_tip/8.2;
        
    end
    attack_s_t = attack_s_t(end);
    
    R_st = - attack_t_design + attack_s_t;
    
    i_s_t = i_opt_rott + R_st / (1+0.5*(Ksh_i*M_R1_t)^3);
    
    i_mt  = i_c_t + (i_s_t - i_c_t)*R_ct/(R_ct+R_st);
    
        %Losses correlations
        
        Km = abs(tand(beta_1_m) - V_2A(end)/V_1A_m * tand(beta_2_m));
        Kt = abs(tand(beta_1_t) - V_2A(end)/V_1A_t * tand(beta_2_t));
        Kh = abs(tand(beta_1_h) - V_2A(end)/V_1A_h * tand(beta_2_h));
        
        Wmax_W1_m = 1.12 + 0.61 * cosd(beta_1_m)^2/sigma_R_m * Km;
        Wmax_W1_t = 1.12 + 0.61 * cosd(beta_1_t)^2/sigma_R_t * Kt;
        Wmax_W1_h = 1.12 + 0.61 * cosd(beta_1_h)^2/sigma_R_h * Kh;
        
        Dm = Wmax_W1_m * W_1_m / W_2_m;
        Dt = Wmax_W1_t * W_1_t / W_2_t;
        Dh = Wmax_W1_h * W_1_h / W_2_h;
        
        Y_2_p_tot_m = 0.004 * ( 1 + 3.1*(Dm-1)^2 + 0.4*(Dm-1)^8 ) * 2 * sigma_R_m / cosd(beta_2_m) * (W_2_m/W_1_m)^2;
        Y_2_p_tot_t = 0.004 * ( 1 + 3.1*(Dt-1)^2 + 0.4*(Dt-1)^8 ) * 2 * sigma_R_t / cosd(beta_2_t) * (W_2_t/W_1_t)^2;
        Y_2_p_tot_h = 0.004 * ( 1 + 3.1*(Dh-1)^2 + 0.4*(Dh-1)^8 ) * 2 * sigma_R_h / cosd(beta_2_h) * (W_2_h/W_1_h)^2;
        
        Y_2_p_tot_m = Y_2_p_tot_m * ( 1 + (i_mm - i_opt_rotm)^2 / R_sm^2 ); 
        Y_2_p_tot_h = Y_2_p_tot_h * ( 1 + (i_mh - i_opt_roth)^2 / R_sh^2 );
        Y_2_p_tot_t = Y_2_p_tot_t * ( 1 + (i_mt - i_opt_rott)^2 / R_st^2 );
        
        % Critical Mach calculation
        
        T_TR1_m = T_1_m + W_1_m^2/2/cp;
        T_TR1_h = T_1_t + W_1_h^2/2/cp;
        T_TR1_t = T_1_t + W_1_t^2/2/cp;
        
        T_cr_1_m = T_TR1_m * 2 / (gamma+1);
        T_cr_1_h = T_TR1_h * 2 / (gamma+1);
        T_cr_1_t = T_TR1_t * 2 / (gamma+1);
        
        W_cr_1_m   = 1 * sqrt(gamma*R_star*T_cr_1_m);
        W_cr_1_h   = 1 * sqrt(gamma*R_star*T_cr_1_h);
        W_cr_1_t   = 1 * sqrt(gamma*R_star*T_cr_1_t);
        
        M_cr_1_m   = M_R1_m * W_cr_1_m / ( W_1_m * Wmax_W1_m );
        M_cr_1_h   = M_R1_h * W_cr_1_h / ( W_1_h * Wmax_W1_h );
        M_cr_1_t   = M_R1_t * W_cr_1_t / ( W_1_t * Wmax_W1_t );
        
        Y_2_p_tot_m = Y_2_p_tot_m + Ksh_i * ((M_R1_m/M_cr_1_m-1)*W_cr_1_m/W_1_m)^2;
        Y_2_p_tot_h = Y_2_p_tot_h + Ksh_i * ((M_R1_h/M_cr_1_h-1)*W_cr_1_h/W_1_h)^2;
        Y_2_p_tot_t = Y_2_p_tot_t + Ksh_i * ((M_R1_t/M_cr_1_t-1)*W_cr_1_t/W_1_t)^2;
        
        %TDN variables
        
        p_TR1_m = p_1_m * ( (1+(gamma-1)/2*(W_1_m^2/(gamma*R_star*T_1_m)))^(gamma/(gamma-1)) );
        p_TR1_t = p_1_t * ( (1+(gamma-1)/2*(W_1_t^2/(gamma*R_star*T_1_t)))^(gamma/(gamma-1)) );
        p_TR1_h = p_1_h * ( (1+(gamma-1)/2*(W_1_h^2/(gamma*R_star*T_1_h)))^(gamma/(gamma-1)) );
        
        p_TR2_m = p_TR1_m - Y_2_p_tot_m * (p_TR1_m - p_1_m);
        p_TR2_t = p_TR1_t - Y_2_p_tot_t * (p_TR1_t - p_1_t);
        p_TR2_h = p_TR1_h - Y_2_p_tot_h * (p_TR1_h - p_1_h);
        
        T_T2_m  = T_T1_m + l_Eu/cp;
        T_T2_t  = T_T1_t + l_Eu/cp;
        T_T2_h  = T_T1_h + l_Eu/cp;
       
        T_2_m = T_T2_m - V_2_m^2 / 2 / cp;
        T_2_t = T_T2_t - V_2_t^2 / 2 / cp;
        T_2_h = T_T2_h - V_2_h^2 / 2 / cp;
        
        p_2_m = p_TR2_m / ( (1+(gamma-1)/2*(W_2_m^2/(gamma*R_star*T_2_m)))^(gamma/(gamma-1)) );
        p_2_t = p_TR2_t / ( (1+(gamma-1)/2*(W_2_t^2/(gamma*R_star*T_2_t)))^(gamma/(gamma-1)) );
        p_2_h = p_TR2_h / ( (1+(gamma-1)/2*(W_2_h^2/(gamma*R_star*T_2_h)))^(gamma/(gamma-1)) );
        
        p_T2_m = p_2_m * ( (1+(gamma-1)/2*(rho_2_m*V_2_m^2/(gamma*p_2_m)))^(gamma/(gamma-1)) );
        p_T2_t = p_2_t * ( (1+(gamma-1)/2*(rho_2_t*V_2_t^2/(gamma*p_2_t)))^(gamma/(gamma-1)) );
        p_T2_h = p_2_h * ( (1+(gamma-1)/2*(rho_2_h*V_2_h^2/(gamma*p_2_h)))^(gamma/(gamma-1)) );
        
        T_2_m_is  = T_1_m * ( p_2_m / p_1_m )^((gamma-1)/gamma);
        T_2_t_is  = T_1_t * ( p_2_t / p_1_t )^((gamma-1)/gamma);
        T_2_h_is  = T_1_h * ( p_2_h / p_1_h )^((gamma-1)/gamma);
        
        T_T2_m_is  = T_T1_m * ( p_T2_m / p_T1_m )^((gamma-1)/gamma);
        T_T2_t_is  = T_T1_t * ( p_T2_t / p_T1_t )^((gamma-1)/gamma);
        T_T2_h_is  = T_T1_h * ( p_T2_h / p_T1_h )^((gamma-1)/gamma);
        
        %Update SS efficiency @ rotor outlet
        
        eta_R_m(end+1) = ( T_2_m_is - T_1_m ) / ( T_2_m - T_1_m );
        eta_R_t(end+1) = ( T_2_t_is - T_1_t ) / ( T_2_t - T_1_t );
        eta_R_h(end+1) = ( T_2_h_is - T_1_h ) / ( T_2_h - T_1_h );
        
        %Update velocity triangle @ rotor outlet
        
            rho_2_m = p_2_m / R_star / T_2_m; 
            rho_2_t = p_2_t / R_star / T_2_t; 
            rho_2_h = p_2_h / R_star / T_2_h; 
            
            V_2A(end+1) = 3 * m / pi / b / D_m / (rho_2_t + rho_2_m + rho_2_h);

    
        end
        
           V_2A_m = V_2A(end);
           V_2A_h = V_2A_h(end);
           V_2A_t = V_2A_t(end);
        
        %Now we optimize the solidity with the Howell correlation and we
        %evaluate again the efficiency
        
           if i==1 && HOW_OPT == 1
               
           Db_Psi = ppval(Dbeta_Psi_curve, abs(beta_2_m));
           Psi_opt = abs(beta_2_m-beta_1_m)/Db_Psi;        
           if Psi_opt<1.329 && Psi_opt>0.749
           x=0.4:0.001:1.6;
           s_over_c_R = mean(x( find(ppval(Psi_curve, x)>Psi_opt-0.001 & ppval(Psi_curve, x)<Psi_opt+0.001)));
           Howell_R='Rotor Optimized';
           elseif Psi_opt>1.329
           s_over_c_R = 0.4;
           Howell_R='Rotor Overloaded';
           else
           s_over_c_R = 1.6;        
           Howell_R='Rotor Underloaded';
           end
           sigma_R_m = 1/s_over_c_R;
           eta_R_t = [2*eta_R_t(end) eta_R_t(end)];
           eta_R_m = [2*eta_R_m(end) eta_R_m(end)];
           eta_R_h = [2*eta_R_h(end) eta_R_h(end)];
           
           end
        
        end
        
        
        story_eta_R_t = eta_R_t;
        story_eta_R_m = eta_R_m; 
        story_eta_R_h = eta_R_h; 
        story_V_2A    = V_2A;

        eta_R_m   = eta_R_m(end);
        eta_R_t   = eta_R_t(end);
        eta_R_h   = eta_R_h(end);
        V_2A      = V_2A(end);
        
        
           

