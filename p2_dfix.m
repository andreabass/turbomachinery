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

        % To initialize the velocity triangle at the outlet of the rotor
        % we refer to an equivalent problem in which the losses are 
        % assumed through the condition eta_R,j = eta_TT for j = m,t,h

        p2_init
        
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
         N_R = floor(pi * D_m / s_R_m);
         
         
           

