
while abs(V_3A(end) - V_3A(end-1))> tol
               V_3A(end-1) = V_3A(end);

        V_3_t = sqrt(V_3T_t^2 + V_3A(end)^2);
        V_3_m = sqrt(V_3T_m^2 + V_3A(end)^2);
        V_3_h = sqrt(V_3T_h^2 + V_3A(end)^2);
        
        alpha_3_t = atand(V_3T_t / V_3A(end));
        alpha_3_m = atand(V_3T_m / V_3A(end));
        alpha_3_h = atand(V_3T_h / V_3A(end));
        
        T_T3_m = T_T2_m;
        T_T3_t = T_T2_t;
        T_T3_h = T_T2_h;
  
        T_3_m = T_T3_m - V_3_m^2 / 2 / cp;
        T_3_t = T_T3_t - V_3_t^2 / 2 / cp;
        T_3_h = T_T3_h - V_3_h^2 / 2 / cp;
        
        T_3_t_is = T_2_t + eta_S_t(end) * (T_3_t - T_2_t);
        T_3_m_is = T_2_m + eta_S_m(end) * (T_3_m - T_2_m);
        T_3_h_is = T_2_h + eta_S_h(end) * (T_3_h - T_2_h);
        
        p_3_t = p_2_t * (T_2_t / T_3_t_is) ^ (gamma/(1-gamma));
        p_3_m = p_2_m * (T_2_m / T_3_m_is) ^ (gamma/(1-gamma));
        p_3_h = p_2_h * (T_2_h / T_3_h_is) ^ (gamma/(1-gamma));
        
        rho_3_t = p_3_t / R_star / T_3_t;
        rho_3_m = p_3_m / R_star / T_3_m;
        rho_3_h = p_3_h / R_star / T_3_h;
        
        V_3A(end+1) = 3 * m / pi / b / D_m / (rho_3_t + rho_3_m + rho_3_h);
end