
% OPTIMUM INCIDENCE AND DEVIATION (ALREADY COMPUTED)

% DESIGN INCIDENCE
i_des_statt = i_S_mt;
i_des_statm = i_S_mm;
i_des_stath = i_S_mh;

% DESIGN DEVIATION ANGLE
ddelta_di_t = ( 1+(sigma_S_t+0.25*sigma_S_t^4)*(abs(alpha_2_t)/53)^2.5 )/exp(3.1*sigma_S_t);
ddelta_di_m = ( 1+(sigma_S_m+0.25*sigma_S_m^4)*(abs(alpha_2_m)/53)^2.5 )/exp(3.1*sigma_S_m);
ddelta_di_h = ( 1+(sigma_S_h+0.25*sigma_S_h^4)*(abs(alpha_2_h)/53)^2.5 )/exp(3.1*sigma_S_h);

dev_des_statt = dev_opt_statt + ddelta_di_t*(i_S_mt - i_opt_statt) + 10*(1-V_3A(end)/V_2A);
dev_des_statm = dev_opt_statm + ddelta_di_m*(i_S_mm - i_opt_statm) + 10*(1-V_3A(end)/V_2A);
dev_des_stath = dev_opt_stath + ddelta_di_h*(i_S_mh - i_opt_stath) + 10*(1-V_3A(end)/V_2A);

% INLET GEOMETRICAL ANGLE
alpha_2_h_geo = alpha_2_h + i_des_stath;
alpha_2_m_geo = alpha_2_m + i_des_statm;
alpha_2_t_geo = alpha_2_t + i_des_statt;

% OUTLET GEOMETRICAL ANGLE
alpha_3_h_geo = alpha_3_h + dev_des_stath; 
alpha_3_m_geo = alpha_3_m + dev_des_statm;
alpha_3_t_geo = alpha_3_t + dev_des_statt;

% CAMBER
teta_hub_stat = abs(alpha_3_h_geo - alpha_2_h_geo); 
teta_mid_stat = abs(alpha_3_m_geo - alpha_2_m_geo); 
teta_tip_stat = abs(alpha_3_t_geo - alpha_2_t_geo); 

% STAGGER ANGLE
gamma_stath =  alpha_3_h_geo - teta_hub_stat/2;
gamma_statm =  alpha_3_m_geo - teta_mid_stat/2;
gamma_statt =  alpha_3_t_geo - teta_tip_stat/2;

% ATTACK ANGLE
attack_S_h_design = gamma_stath - alpha_2_h;
attack_S_m_design = gamma_statm - alpha_2_m;
attack_S_t_design = gamma_statt - alpha_2_t;