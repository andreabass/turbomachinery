        mass_IGV   =  (rho_0(end) * V_0A * pi * b * D_m) - (mean(rho_1) * V_1A * pi * b * D_m);
        energy_IGV = [ T_T0 - (T_1_t + V_1_t^2/2/cp); T_T0 - (T_1_m + V_1_m^2/2/cp); T_T0 - (T_1_h + V_1_h^2/2/cp)];
       
        mass_R   = ( mean(rho_1) * V_1A * pi * b *D_m ) -  ( mean(rho_2) * V_2A * pi * b * D_m );
        energy_R = [ l_Eu - cp*(T_T2_m - T_T1_m) ; l_Eu - cp*(T_T2_t - T_T1_t); l_Eu - cp*(T_T2_h - T_T1_h)];
       
        
        mass_S   = ( mean(rho_2) * V_2A * pi * b *D_m ) - ( mean(rho_3) * V_3A(end) * pi * b * D_m );
        energy_S = [ T_T2_t - (T_3_t + V_3_t^2/2/cp); T_T2_m - (T_3_m + V_3_m^2/2/cp); T_T2_h - (T_3_h + V_3_h^2/2/cp)];
        
        mass = [mass_IGV, mass_R, mass_S];
        energy = [energy_IGV, energy_R, energy_S ];
        
        disp('__________________________________')
        disp(' ')
        disp('             CHECKS')
        disp('__________________________________')
        disp(' ')
        
        if nnz(energy<1e-6)==numel(energy) && nnz(mass<1e-6)==numel(mass)
            disp('[ OK ] Mass and energy balances')
        else
            disp('[ ERROR ] Mass and energy balances')
        end
        
        if beta_TT*0.99999 <= (p_T3_av / p_T0) <= beta_TT*1.00001
    disp('[ OK ] Compression ratio')
        else
            
    disp('[ ERROR ] Compression ratio')
        end

if m*0.99999 <= mean(rho_1) * V_1A * pi * b *D_m <= m*1.00001
    disp('[ OK ] Discharged mass flow rate')
else
        disp('[ ERROR ] Discharged mass flow rate')
end

disp(' ')




            
       
    
  