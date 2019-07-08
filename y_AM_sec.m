function Y = y_AM_sec(alpha_out,alpha_in,c,b,s_c)
        
        alpha_av_t_01 = atand((tand(alpha_in)+tand(alpha_out))/2);
        cL = 2 * s_c * (abs(tand(alpha_out)-tand(alpha_in)))*cosd(alpha_av_t_01);
    
        Y = c / b *(0.0334 * cosd(alpha_out)/cosd(alpha_in)) * (cL / s_c)^2 * ((cosd(alpha_out))^2) / (cosd(alpha_av_t_01))^3;
    
end

