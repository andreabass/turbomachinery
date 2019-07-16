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

        % We initialize the total pressures at the stator outlet at the 
        % value of the average total pressure required by the client
        
        % We initialize the outlet velocity field assuming constant axial
        % velocity
        
        p_T3_m = [beta_TT(end)*p_T0 beta_TT(end)*p_T0*(1+tol)];
        p_T3_t = [beta_TT(end)*p_T0 beta_TT(end)*p_T0*(1+tol)];
        p_T3_h = [beta_TT(end)*p_T0 beta_TT(end)*p_T0*(1+tol)];   
        
        V_3A   = V_2A;
        
        while abs(p_T3_m(end) - p_T3_m(end-1))>tol || abs(p_T3_t(end) - p_T3_t(end-1))>tol || abs(p_T3_h(end) - p_T3_h(end-1))>tol
            p_T3_m(end-1) = p_T3_m(end);
            p_T3_t(end-1) = p_T3_t(end);
            p_T3_h(end-1) = p_T3_h(end);
            
        V_3_m = sqrt( V_3T_m^2 + V_3A^2 );
        V_3_t = sqrt( V_3T_t^2 + V_3A^2 );
        V_3_h = sqrt( V_3T_h^2 + V_3A^2 );
        
        alpha_3_m = atand(V_3T_m / V_3A);
        alpha_3_t = atand(V_3T_t / V_3A);
        alpha_3_h = atand(V_3T_h / V_3A);
        
        % Energy balance
        
        T_T3_m = T_T2_m;
        T_T3_t = T_T2_t;
        T_T3_h = T_T2_h;
        
        % Static temperature field
        
        T_3_m = T_T3_m - V_3_m^2 / 2 / cp;
        T_3_t = T_T3_t - V_3_t^2 / 2 / cp;
        T_3_h = T_T3_h - V_3_h^2 / 2 / cp;
        
        p_3_m = p_T3_m(end) / ( 1 + V_3_m^2 / 2 / R_star / T_3_m );
        p_3_t = p_T3_t(end) / ( 1 + V_3_t^2 / 2 / R_star / T_3_t );
        p_3_h = p_T3_h(end) / ( 1 + V_3_h^2 / 2 / R_star / T_3_h );
        
        rho_3_m = p_3_m / R_star / T_3_m ;
        rho_3_t = p_3_t / R_star / T_3_t ;
        rho_3_h = p_3_h / R_star / T_3_h ;
        
        %   HOWELL CORRELATION 
        % These assumptions are used to obtain the two correction
        % coefficients of Howell correlation (Psi, Phi) equal to one. In
        % principle is necessary to check also that the deflection Dbeta is
        % sufficiently close to the optimal value associated to beta2
        
        %    Psi_opt = Db_Psi/(alpha_2_m-alpha_1_m);
        %    x=0.4:0.001:1.6;
        %    s_over_c_S = mean(x( find(ppval(Psi_curve, x)>Psi_opt-0.001 & ppval(Psi_curve, x)<Psi_opt+0.001)));       
        % sigma_S_m = 1/s_over_c_S;
        
        sigma_S_m = 1;
        
        %Chord calculation based on Howell value for Reynolds number
        
            c_How_S_m = mu*Re_How/rho_3_m/V_3_m;
            c_How_S_t = mu*Re_How/rho_3_t/V_3_t;
            c_How_S_h = mu*Re_How/rho_3_h/V_3_h;
         c_S_m = 1.1*max(max(c_How_S_m,c_How_S_t), c_How_S_h);
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
        
        p_T3_m(end+1) = p_T2_m - Y_3_p_tot_m * (p_T2_m - p_2_m);
        p_T3_t(end+1) = p_T2_t - Y_3_p_tot_t * (p_T2_t - p_2_t);
        p_T3_h(end+1) = p_T2_h - Y_3_p_tot_h * (p_T2_h - p_2_h);      

        end
        
        
        
        
           

