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
        Re_How = 3e5;
        sigma_R_m = 1;
        
        %Chord calculation based on Howell value for Reynolds number
         

