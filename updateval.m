
    p_T3_av = 0.3*p_T3_h + 0.4*p_T3_m + 0.3*p_T3_t;
    T_T3_av = 0.3*T_T3_h + 0.4*T_T3_m + 0.3*T_T3_t;
    
    eta_TT(end+1) = cp * T_T0 * ( (p_T3_av/p_T0) ^ ((gamma - 1)/gamma) - 1 ) / (cp *  (T_T3_av - T_T0) );
    
    
    

    