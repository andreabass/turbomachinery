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
        
        eta_R_t = [eta_TT(end) eta_TT(end) + 2*tol];
        eta_R_m = [eta_TT(end) eta_TT(end) + 2*tol];
        eta_R_h = [eta_TT(end) eta_TT(end) + 2*tol];    
        
        % Constant axial velocity as initial value
        
        V_2A = V_1A_m;
        V_2A = [V_2A V_2A + 2*tol];
        Howell_correlation
        opt=1;
        
        for i=1:2
        while abs(eta_R_t(end) - eta_R_t(end-1))>tol || abs(eta_R_m(end) - eta_R_m(end-1))>tol || abs(eta_R_h(end) - eta_R_h(end-1))>tol
            eta_R_t(end-1) = eta_R_t(end);
            eta_R_m(end-1) = eta_R_m(end);
            eta_R_h(end-1) = eta_R_h(end);
        
            p2_rotvelout
        
        %   HOWELL CORRELATION %
        %These assumptions are used to obtain the two correction
        %coefficients of Howell correlation (Psi, Phi) equal to one. In
        %principle is necessary to check also that the deflection Dbeta is
        %sufficiently close to the optimal value associated to beta2
            
        Re_How = 3e5;
       
        if i==1
            sigma_R_m = 1; 
        end
        
        
        %Chord calculation based on Howell value for Reynolds number
        
            c_How_m = mu*Re_How/rho_2_m/W_2_m;
            c_How_t = mu*Re_How/rho_2_t/W_2_t;
            c_How_h = mu*Re_How/rho_2_h/W_2_h;
         c_R_m = 1.1*max(max(c_How_m,c_How_t), c_How_h);
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
        
        %Losses correlations
        
        Dm = (W_1_m-W_2_m)/W_1_t+abs((W_1T_m-W_2T_m)/(2*W_1_m*sigma_R_m));
        Dt = (W_1_t-W_2_t)/W_1_t+abs((W_1T_t-W_2T_t)/(2*W_1_t*sigma_R_t));
        Dh = (W_1_h-W_2_h)/W_1_h+abs((W_1T_h-W_2T_h)/(2*W_1_h*sigma_R_h));
        
        Y_2_p_tot_m = 0.0035*(1+3.5*Dm+37*(Dm)^4)*2*sigma_R_m/cosd(beta_2_m);
        Y_2_p_tot_t = 0.0035*(1+3.5*Dt+37*(Dt)^4)*2*sigma_R_t/cosd(beta_2_t);
        Y_2_p_tot_h = 0.0035*(1+3.5*Dh+37*(Dh)^4)*2*sigma_R_h/cosd(beta_2_h);
        
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
        
        %Now we optimize the solidity with the Howell correlation and we
        %evaluate again the efficiency
        
           if i==1
               
           Db_Psi = ppval(Dbeta_Psi_curve, abs(beta_2_m));
           Psi_opt = Db_Psi/(beta_2_m-beta_1_m);        
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
           eta_R_t = [eta_R_t(end) eta_R_t(end) + 2*tol];
           eta_R_m = [eta_R_m(end) eta_R_m(end) + 2*tol];
           eta_R_h = [eta_R_h(end) eta_R_h(end) + 2*tol];
           
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
        
        
           

