%3rd subproblem (VGV outlet / Rotor inlet thermodynamics and losses

T_1_t = T_T1_t - (V_1_t^2) / (2*cp);
T_1_m = T_T1_m - (V_1_m^2) / (2*cp);
T_1_h = T_T1_h - (V_1_h^2) / (2*cp);

t_over_s_t = 0.02;
t_over_s_m = 0.02;
t_over_s_h = 0.02;

%Losses calculation
alpha_2prime_t = 90 - alpha_1_t;
alpha_2prime_m = 90 - alpha_1_m;
alpha_2prime_h = 90 - alpha_1_h;

%Mid
if alpha_2prime_m > 30
    s_over_c_min = 0.614 + alpha_2prime_m / 130;
else
    s_over_c_min = 0.46 + alpha_2prime_m / 77;
end

if alpha_2prime_m > 27
    A = 0.025 + (27 - alpha_2prime_m) / 3085;
else
    A = 0.025 + (27 - alpha_2prime_m) / 530;
end

Y_p_1_in = A;
Re_ref = 200000;
alpha_av_t_01 = atand((tand(alpha_0_m)+tand(alpha_1_m))/2);
cL = 2 * s_over_c_min * (abs(tand(alpha_1_m)-tand(alpha_0_m)))*cosd(alpha_av_t_01);
c_IGV = 0.04;
c_IGV_t = c_IGV;
c_IGV_m = c_IGV;
c_IGV_h = c_IGV;
Re_1_m = V_1_m * rho_1_m(end) * c_IGV_m / mu;
Y_p_1_Re = Y_p_1_in * (Re_ref / Re_1_m)^0.2;
Y_1_sec = c_IGV_m / b *(0.0334 * cosd(alpha_1_m)/cosd(alpha_0_m)) * (cL / s_over_c_min)^2 * ((cosd(alpha_1_m))^2) / (cosd(alpha_av_t_01))^3;
Y_1_p_tot = Y_p_1_Re + Y_1_sec;

p_T1_m = p_T0_m - Y_1_p_tot * (p_T0_m - p_0_m);
p_1_m = p_T1_m / (1 + (V_1_m^2)/(2 * R_star * T_1_m));

rho_1_m(end) = p_1_m / R_star / T_1_m;

s_IGV_m = s_over_c_min * c_IGV;
N_bl_IGV = pi * D_m / s_IGV_m;
s_IGV_t = pi * D_t / N_bl_IGV;
s_IGV_h = pi * D_h / N_bl_IGV;

s_over_c_t = s_IGV_t / c_IGV;
s_over_c_h = s_IGV_h / c_IGV;

if alpha_2prime_t > 27
    A_t = 0.025 + (27 - alpha_2prime_t) / 3085;
else
    A_t = 0.025 + (27 - alpha_2prime_t) / 530;
end

if alpha_2prime_h > 27
    A_h = 0.025 + (27 - alpha_2prime_h) / 3085;
else
    A_h = 0.025 + (27 - alpha_2prime_h) / 530;
end
%Tip
Y_p_1_in_t = A_t;
alpha_av_t_01_t = atand((tand(alpha_0_t)+tand(alpha_1_t))/2);
cL_t = 2 * s_over_c_t * (abs(tand(alpha_1_t)-tand(alpha_0_t)))*cosd(alpha_av_t_01_t);
Re_1_t = V_1_t * rho_1_t(end) * c_IGV_t / mu;
Y_p_1_Re_t = Y_p_1_in_t * (Re_ref / Re_1_t)^0.2;
Y_1_sec_t = c_IGV_t / b *(0.0334 * cosd(alpha_1_t)/cosd(alpha_0_t)) * (cL_t / s_over_c_t)^2 * ((cosd(alpha_1_t))^2) / (cosd(alpha_av_t_01_t))^3;
Y_1_p_tot_t = Y_p_1_Re_t + Y_1_sec_t;
p_T1_t = p_T0_t - Y_1_p_tot_t * (p_T0_t - p_0_t);
p_1_t = p_T1_t / (1 + (V_1_t^2)/(2 * R_star * T_1_t));

rho_1_t(end) = p_1_t / R_star / T_1_t;

%Hub
Y_p_1_in_h = A_h;
alpha_av_t_01_h = atand((tand(alpha_0_h)+tand(alpha_1_h))/2);
cL_h = 2 * s_over_c_h * (abs(tand(alpha_1_h)-tand(alpha_0_h)))*cosd(alpha_av_t_01_h);
Re_1_h = V_1_h * rho_1_h(end) * c_IGV_h / mu;
Y_p_1_Re_h = Y_p_1_in_h * (Re_ref / Re_1_h)^0.2;
Y_1_sec_h = c_IGV_h / b *(0.0334 * cosd(alpha_1_h)/cosd(alpha_0_h)) * (cL / s_over_c_h)^2 * ((cosd(alpha_1_h))^2) / (cosd(alpha_av_t_01_h))^3;
Y_1_p_tot_h = Y_p_1_Re_h + Y_1_sec_h;
p_T1_h = p_T0_h - Y_1_p_tot_h * (p_T0_h - p_0_h);
p_1_h = p_T1_h / (1 + (V_1_h^2)/(2 * R_star * T_1_h));

rho_1_h(end) = p_1_h / R_star / T_1_h;