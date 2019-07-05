         %% ROTOR  OUTLET %%
          %%%%%%%%%%%%%%%
           %%%%%%%%%%%%%
            %%%%%%%%%%%
             %%%%%%%%%
              %%%%%%%
               %%%%%
                %%%
                 %
        
                 
        %%%%%% [INITIALIZATION] %%%%%%%
        %Initialization value for V_2A_m
        V_2A_m = V_1A_m;
        V_2A_t = V_2A_m;
        V_2A_h = V_2A_m;
        
        % Let's define initialization values for inlet axial velocity and
        % rotor efficiency defined as the ratio between isentropic static
        % enthalpy change and real static enthalpy change. Thanks to these
        % values we'll be able to fully determine the VT and TDN quantities
        % at the rotor outlet that will allow us to calculate losses and
        % therefore get to adjusted values of the TDN quantities. Such
        % updated values will be finally used to recalculate the efficiency
        % of the rotor.In the meanwhile for every cycle we'll impose the
        % mass balance to update the axial velocity.
        
        V_2A = [V_2A_m V_2A_m + 2*tol];
        eta_R_t = [eta_TT_m - 2*tol eta_TT_m];
        eta_R_m = [eta_TT_m - 2*tol eta_TT_m];
        eta_R_h = [eta_TT_m - 2*tol eta_TT_m];
        
        l_Eu = deltaHis_TT / eta_TT_m;
%        while abs(eta_R_t(end) - eta_R_t(end-1))> tol || abs(eta_R_m(end) - eta_R_m(end-1))> tol || abs(eta_R_h(end) - eta_R_h(end-1))> tol
        
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
         
         
           

