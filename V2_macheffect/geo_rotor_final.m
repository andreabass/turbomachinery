
% OPTIMUM INCIDENCE AND DEVIATION (ALREADY COMPUTED)

% DESIGN INCIDENCE
i_des_rott = i_mt;
i_des_rotm = i_mm;
i_des_roth = i_mh;

% DESIGN DEVIATION ANGLE
ddelta_di_t = ( 1+(sigma_R_t+0.25*sigma_R_t^4)*(abs(beta_1_t)/53)^2.5 )/exp(3.1*sigma_R_t);
ddelta_di_m = ( 1+(sigma_R_m+0.25*sigma_R_m^4)*(abs(beta_1_m)/53)^2.5 )/exp(3.1*sigma_R_m);
ddelta_di_h = ( 1+(sigma_R_h+0.25*sigma_R_h^4)*(abs(beta_1_h)/53)^2.5 )/exp(3.1*sigma_R_h);
    
dev_des_rott = dev_opt_rott + ddelta_di_t*(i_mt - i_opt_rott) + 10*(1-V_2A(end)/V_1A);
dev_des_rotm = dev_opt_rotm + ddelta_di_m*(i_mm - i_opt_rotm) + 10*(1-V_2A(end)/V_1A);
dev_des_roth = dev_opt_roth + ddelta_di_h*(i_mh - i_opt_roth) + 10*(1-V_2A(end)/V_1A);

% INLET GEOMETRICAL ANGLE
beta_1_h_geo = beta_1_h + i_des_roth;
beta_1_m_geo = beta_1_m + i_des_rotm;
beta_1_t_geo = beta_1_t + i_des_rott;

% OUTLET GEOMETRICAL ANGLE
beta_2_h_geo = beta_2_h + dev_des_roth; 
beta_2_m_geo = beta_2_m + dev_des_rotm;
beta_2_t_geo = beta_2_t + dev_des_rott;

% CAMBER
teta_hub = abs(beta_2_h_geo - beta_1_h_geo); 
teta_mid = abs(beta_2_m_geo - beta_1_m_geo); 
teta_tip = abs(beta_2_t_geo - beta_1_t_geo); 

% STAGGER ANGLE
gamma_roth =  beta_2_h_geo - teta_hub/2;
gamma_rotm =  beta_2_m_geo - teta_mid/2;
gamma_rott =  beta_2_t_geo - teta_tip/2;

% ATTACK ANGLE
attack_h_design = gamma_roth - beta_1_h;
attack_m_design = gamma_rotm - beta_1_m;
attack_t_design = gamma_rott - beta_1_t;
