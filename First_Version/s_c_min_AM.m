function s_over_c_min = s_c_min_AM(alpha_out)

        alpha_2prime = 90 - alpha_out;
        if alpha_2prime > 30
        s_over_c_min = 0.614 + alpha_2prime / 130;
        else
   s_over_c_min = 0.46 + alpha_2prime / 77;
        end

end

