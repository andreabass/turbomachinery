function Y = y_AM_inc_min(alpha_out)

        alpha_2prime = 90 - alpha_out;

        if alpha_2prime > 27
        A = 0.025 + (27 - alpha_2prime) / 3085;
        else
        A = 0.025 + (27 - alpha_2prime) / 530;
        end

        Y = A;

end

