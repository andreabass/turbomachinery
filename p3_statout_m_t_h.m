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
        
        V_3T_m = 0;
        
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
        
        while abs(eta_S_t(end) - eta_S_t(end-1))>tol || abs(eta_S_m(end) - eta_S_m(end-1))>tol || abs(eta_S_h(end) - eta_S_h(end-1))>tol
            eta_S_t(end-1) = eta_S_t(end);
            eta_S_m(end-1) = eta_S_m(end);
            eta_S_h(end-1) = eta_S_h(end);
        
            p3_statvelout
        
        V_3A_t = V_3A;
        V_3A_h = V_3A;
        
        %   HOWELL CORRELATION %
        %These assumptions are used to obtain the two correction
        %coefficients of Howell correlation (Psi, Phi) equal to one. In
        %principle is necessary to check also that the deflection Dbeta is
        %sufficiently close to the optimal value associated to beta2
            
        Re_How = 3e5;
       
        %    Howell_correlation
        %    Db_Psi = ppval(Dbeta_Psi_curve, abs(alpha_3_m));
        %    Psi_opt = Db_Psi/(alpha_3_m-beta_1_m);
        %    x=0.4:0.001:1.6;
        %    s_over_c = mean(x( find(ppval(Psi_curve, x)>Psi_opt-0.001 & ppval(Psi_curve, x)<Psi_opt+0.001)));       
        % sigma_R_m = 1/s_over_c;
        
        sigma_S_m = 1; 
        
        %Chord calculation based on Howell value for Reynolds number
        
            c_How_m_S = mu*Re_How/rho_3_m/V_3_m;
            c_How_t_S = mu*Re_How/rho_3_t/V_3_t;
            c_How_h_S = mu*Re_How/rho_3_h/V_3_h;
         c_S_m = 1.1*max(max(c_How_m_S,c_How_t_S), c_How_h_S);
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
        
        %Losses correlations
        
        DmS = (V_2_m-V_3_m)/V_2_t+abs((V_2T_m-V_3T_m)/(2*V_2_m*sigma_S_m));
        DtS = (V_2_t-V_3_t)/V_2_t+abs((V_2T_t-V_3T_t)/(2*V_2_t*sigma_S_t));
        DhS = (V_2_h-V_3_h)/V_2_h+abs((V_2T_h-V_3T_h)/(2*V_2_h*sigma_S_h));
        
        Y_3_p_tot_m = 0.0035*(1+3.5*DmS+37*(DmS)^4)*2*sigma_S_m/cosd(alpha_3_m);
        Y_3_p_tot_t = 0.0035*(1+3.5*DtS+37*(DtS)^4)*2*sigma_S_t/cosd(alpha_3_t);
        Y_3_p_tot_h = 0.0035*(1+3.5*DhS+37*(DhS)^4)*2*sigma_S_h/cosd(alpha_3_h);
        
        %TDN variables
        
        p_T3_m = p_T2_m - Y_3_p_tot_m * (p_T2_m - p_2_m);
        p_T3_t = p_T2_t - Y_3_p_tot_t * (p_T2_t - p_2_t);
        p_T3_h = p_T2_h - Y_3_p_tot_h * (p_T2_h - p_2_h);  
        
        p_3_m   = p_T3_m - rho_3_m * V_3_m^2 / 2;
        p_3_t   = p_T3_t - rho_3_t * V_3_t^2 / 2;
        p_3_h   = p_T3_h - rho_3_h * V_3_h^2 / 2;
        
        T_3_m_is  = T_2_m * ( p_3_m / p_2_m )^((gamma-1)/gamma);
        T_3_t_is  = T_2_t * ( p_3_t / p_2_t )^((gamma-1)/gamma);
        T_3_h_is  = T_2_h * ( p_3_h / p_2_h )^((gamma-1)/gamma);
        
        T_T3_m_is  = T_T2_m * ( p_T3_m / p_T2_m )^((gamma-1)/gamma);
        T_T3_t_is  = T_T2_t * ( p_T3_t / p_T2_t )^((gamma-1)/gamma);
        T_T3_h_is  = T_T2_h * ( p_T3_h / p_T2_h )^((gamma-1)/gamma);
        
        T_3_m = T_T3_m - V_3_m^2 / 2 / cp;
        T_3_t = T_T3_t - V_3_t^2 / 2 / cp;
        T_3_h = T_T3_h - V_3_h^2 / 2 / cp;
        
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
        
        story_eta_S_t = eta_S_t;
        story_eta_S_m = eta_S_m; 
        story_eta_S_h = eta_S_h; 
        story_V_3A    = V_3A;

        eta_S_m   = eta_S_m(end);
        eta_S_t   = eta_S_t(end);
        eta_S_h   = eta_S_h(end);
        V_3      = V_3A(end);
        
        
        
        
           

