        mass_IGV   =  (rho_0(end) * V_0A * pi * b * D_m) - (mean(rho_1) * V_1A * pi * b * D_m);
        energy_IGV = [ T_T0 - (T_1_t + V_1_t^2/2/cp); T_T0 - (T_1_m + V_1_m^2/2/cp); T_T0 - (T_1_h + V_1_h^2/2/cp)];
       
        mass_R   = ( mean(rho_1) * V_1A * pi * b *D_m ) -  ( mean(rho_2) * V_2A * pi * b * D_m );
        energy_R = [ l_Eu - cp*(T_T2_m - T_T1_m) ; l_Eu - cp*(T_T2_t - T_T1_t); l_Eu - cp*(T_T2_h - T_T1_h)];
       
        
        mass_S   = ( mean(rho_2) * V_2A * pi * b *D_m ) - ( mean(rho_3) * V_3A(end) * pi * b * D_m );
        energy_S = [ T_T2_t - (T_3_t + V_3_t^2/2/cp); T_T2_m - (T_3_m + V_3_m^2/2/cp); T_T2_h - (T_3_h + V_3_h^2/2/cp)];
        
        mass = [mass_IGV, mass_R, mass_S];
        energy = [energy_IGV, energy_R, energy_S ];
        
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

if all(M_1 < 1)
    disp('Subsonic flow at the IGV outlet (section 1)')
else
        disp('Supersonic flow at the IGV outlet (section 1)')
        
        if all(M_1 < 1.4)
            disp('    Mach number at IGV outlet lower than 1.4')
        end
end

if all(MR_2 < 1)
    disp('Subsonic flow at the rotor outlet (section 2)')
else
        disp('Supersonic flow at the rotor outlet (section 2)')
                if all(MR_2 < 1.4)
            disp('    Relative Mach number at rotor outlet lower than 1.4')
        end
       
end

if all(M_3 < 1)
    disp('Subsonic flow at the stator outlet (section 3)')
else
        disp('Supersonic flow at the stator outlet (section 3)')
                if all(M_3 < 1.4)
            disp('    Mach number at stator outlet lower than 1.4')
        end
end


            
       
    
  