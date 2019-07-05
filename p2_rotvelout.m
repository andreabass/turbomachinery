       
        V_2T_m = V_1T_m + l_Eu / U_m;
        
        V_2T_t = V_2T_m * D_m / D_t;
        V_2T_h = V_2T_m * D_m / D_h(end);
        
        W_2T_m = V_2T_m - U_m;
        W_2T_t = V_2T_t - U_t;
        W_2T_h = V_2T_h - U_h;

while abs(V_2A(end) - V_2A(end-1))> tol
               V_2A(end-1) = V_2A(end);
        
        W_2A_m = V_2A(end);
        W_2A_t = V_2A(end);
        W_2A_h = V_2A(end);
        
               W_2A = W_2A_m;
        
        V_2_t = sqrt(V_2T_t^2 + V_2A(end)^2);
        V_2_m = sqrt(V_2T_m^2 + V_2A(end)^2);
        V_2_h = sqrt(V_2T_h^2 + V_2A(end)^2);
        
        W_2_t = sqrt(W_2T_t^2 + W_2A^2);
        W_2_m = sqrt(W_2T_m^2 + W_2A^2);
        W_2_h = sqrt(W_2T_h^2 + W_2A^2);
        
        alpha_2_t = atand(V_2T_t / V_2A(end));
        alpha_2_m = atand(V_2T_m / V_2A(end));
        alpha_2_h = atand(V_2T_h / V_2A(end));
        
        beta_2_t = atand(W_2T_t / W_2A);
        beta_2_m = atand(W_2T_m / W_2A);
        beta_2_h = atand(W_2T_h / W_2A);
        
        Xi_t = (W_1_t^2 - W_2_t^2) / 2 / l_Eu;
        Xi_m = (W_1_m^2 - W_2_m^2) / 2 / l_Eu;
        Xi_h = (W_1_h^2 - W_2_h^2) / 2 / l_Eu;
        
        DeltaH_R_t = Xi_t * l_Eu;
        DeltaH_R_m = Xi_m * l_Eu;
        DeltaH_R_h = Xi_h * l_Eu;
        
        T_2_t = T_1_t + DeltaH_R_t / cp;
        T_2_m = T_1_m + DeltaH_R_m / cp;
        T_2_h = T_1_h + DeltaH_R_h / cp;
        
        T_2_t_is = T_1_t + eta_R_t(end) * (T_2_t - T_1_t);
        T_2_m_is = T_1_m + eta_R_m(end) * (T_2_m - T_1_m);
        T_2_h_is = T_1_h + eta_R_h(end) * (T_2_h - T_1_h);
        
        p_2_t = p_1_t * (T_1_t / T_2_t_is) ^ (gamma/(1-gamma));
        p_2_m = p_1_m * (T_1_m / T_2_m_is) ^ (gamma/(1-gamma));
        p_2_h = p_1_h * (T_1_h / T_2_h_is) ^ (gamma/(1-gamma));
        
        rho_2_t = p_2_t / R_star / T_2_t;
        rho_2_m = p_2_m / R_star / T_2_m;
        rho_2_h = p_2_h / R_star / T_2_h;
        
        V_2A(end+1) = 3 * m / pi / b / D_m / (rho_2_t + rho_2_m + rho_2_h);
end