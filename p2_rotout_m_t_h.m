         %% ROTOR  OUTLET %%
          %%%%%%%%%%%%%%%
           %%%%%%%%%%%%%
            %%%%%%%%%%%
             %%%%%%%%%
              %%%%%%%
               %%%%%
                %%%
                 %
            
        %%% INITIALIZATION %%%

        % We initialize the losses in the rotor by assuming each static
        % efficiency (mid, tip. hub) equal to the total-to-total efficiency
        % of the stage
        
        eta_R_t = eta_TT_m;
        eta_R_m = eta_TT_m;
        eta_R_h = eta_TT_m;
        
        eta_R_t = [eta_R_t eta_R_t + 2*tol];
        eta_R_m = [eta_R_m eta_R_m + 2*tol];
        eta_R_h = [eta_R_h eta_R_h + 2*tol];
        
        % while abs(eta_R_t(end) - eta_R_t(end-1))>tol || abs(eta_R_m(end) - eta_R_m(end-1))>tol || abs(eta_R_h(end) - eta_R_h(end-1))>tol
        % eta_R_t(end-1) = eta_R_t(end);
        % eta_R_m(end-1) = eta_R_m(end);
        % eta_R_h(end-1) = eta_R_h(end);
        
        p2_rotvelout
        
        %   HOWELL CORRELATION %
        %These assumptions are used to obtain the two correction
        %coefficients of Howell correlation (Psi, Phi) equal to one. In
        %principle is necessary to check also that the deflection Dbeta is
        %sufficiently close to the optimal value associated to beta2
        Re_How = 3e5;
        sigma_R_m = 1;
        
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
        s_R_h = pi * D_h / N_R;
        sigma_R_h = c_R_h / s_R_h;
        
        %
         
           

