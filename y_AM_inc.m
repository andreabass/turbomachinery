function Y = y_AM_inc(alpha_out,s_c)

        alpha_2prime_t = 90 - alpha_out;
        
        C_t = 0.08*((alpha_2prime_t/30)^2-1);
        n_AM_t = 1+alpha_2prime_t/30;
    
        if alpha_2prime_t > 27
        A_t = 0.025 + (27 - alpha_2prime_t) / 3085;
        else
        A_t = 0.025 + (27 - alpha_2prime_t) / 530;
        end
        
        if alpha_2prime_t > 30
        s_over_c_min_t = 0.614 + alpha_2prime_t / 130;
        B_t = 0; % Not available on Aungier (keep Y_p_1_in_t saturated at A)
        X_AM_t = s_c - s_over_c_min_t;
        Y = A_t + B_t * (abs(X_AM_t))^n_AM_t;
        else
        s_over_c_min_t = 0.46 + alpha_2prime_t / 77;
        B_t = 0.1583-alpha_2prime_t/1640;
        X_AM_t = s_over_c_t - s_over_c_min_t;
        Y = A_t+B_t*X_AM_t^2+C_t*X_AM_t^3;
        end

end

