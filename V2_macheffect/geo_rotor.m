
Dbeta_12_t = abs(beta_1_t-beta_2_t);
Dbeta_12_m = abs(beta_1_m-beta_2_m);
Dbeta_12_h = abs(beta_1_h-beta_2_h);

% Incidence: i_opt_rot = Ksh*Kth*i0_10+n_LIE*teta
Ksh_i = 0.7; % DCA

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma_R_t^3/160;
p_mid = 0.914+sigma_R_m^3/160;
p_hub = 0.914+sigma_R_h^3/160;

i0_10tip = abs(beta_1_t)^p_tip/(5+46*exp(-2.3*sigma_R_t))-0.1*sigma_R_t^3*exp((abs(beta_1_t)-70)/4);
i0_10mid = abs(beta_1_m)^p_mid/(5+46*exp(-2.3*sigma_R_m))-0.1*sigma_R_m^3*exp((abs(beta_1_m)-70)/4);
i0_10hub = abs(beta_1_h)^p_hub/(5+46*exp(-2.3*sigma_R_h))-0.1*sigma_R_h^3*exp((abs(beta_1_h)-70)/4);

n_LIEt = 0.025*sigma_R_t-0.06-(abs(beta_1_t)/90)^(1+1.2*sigma_R_t)/(1.5+0.43*sigma_R_t);
n_LIEm = 0.025*sigma_R_m-0.06-(abs(beta_1_m)/90)^(1+1.2*sigma_R_m)/(1.5+0.43*sigma_R_m);
n_LIEh = 0.025*sigma_R_h-0.06-(abs(beta_1_h)/90)^(1+1.2*sigma_R_h)/(1.5+0.43*sigma_R_h);

i_0_rott = Ksh_i*Kth*i0_10tip;
i_0_rotm = Ksh_i*Kth*i0_10mid;
i_0_roth = Ksh_i*Kth*i0_10hub;

% Deviation: dev_opt_rot = Ksh*Kth*dev0_10+m_LIE*teta

Ksh_d = 0.75;

x_LIEt = abs(beta_1_t)/100;
x_LIEm = abs(beta_1_m)/100;
x_LIEh = abs(beta_1_h)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;
m10_m = 0.17-0.0333*x_LIEm+0.333*x_LIEm^2;
m10_h = 0.17-0.0333*x_LIEh+0.333*x_LIEh^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;
b_m = 0.9625-0.17*x_LIEm-0.85*x_LIEm^3;
b_h = 0.9625-0.17*x_LIEh-0.85*x_LIEh^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma_R_t)^b_t;
m_m = m10_m/(sigma_R_m)^b_m;
m_h = m10_h/(sigma_R_h)^b_h;

dev0_10t = 0.01*sigma_R_t*abs(beta_1_t)+(0.74*sigma_R_t^1.9+3*sigma_R_t)*(abs(beta_1_t)/90)^(1.67+1.09*sigma_R_t);
dev0_10m = 0.01*sigma_R_m*abs(beta_1_m)+(0.74*sigma_R_m^1.9+3*sigma_R_m)*(abs(beta_1_m)/90)^(1.67+1.09*sigma_R_m);
dev0_10h = 0.01*sigma_R_h*abs(beta_1_h)+(0.74*sigma_R_h^1.9+3*sigma_R_h)*(abs(beta_1_h)/90)^(1.67+1.09*sigma_R_h);

dev_0_rott = Ksh_d*Kth_dev*dev0_10t;
dev_0_rotm = Ksh_d*Kth_dev*dev0_10m;
dev_0_roth = Ksh_d*Kth_dev*dev0_10h;

% CAMBER ANGLE
teta_tip = (Dbeta_12_t- i_0_rott + dev_0_rott)/(1-m_t+n_LIEt);
teta_mid = (Dbeta_12_m- i_0_rotm + dev_0_rotm)/(1-m_m+n_LIEm);
teta_hub = (Dbeta_12_h- i_0_roth + dev_0_roth)/(1-m_h+n_LIEh);

% INCIDENCE ANGLE
i_opt_rott = Ksh_i*Kth*i0_10tip+n_LIEt*teta_tip;
i_opt_rotm = Ksh_i*Kth*i0_10mid+n_LIEm*teta_mid;
i_opt_roth = Ksh_i*Kth*i0_10hub+n_LIEh*teta_hub;

% DEVIATION ANGLE
dev_opt_rott = Ksh_d*Kth_dev*dev0_10t+m_t*teta_tip;
dev_opt_rotm = Ksh_d*Kth_dev*dev0_10m+m_m*teta_mid;
dev_opt_roth = Ksh_d*Kth_dev*dev0_10h+m_h*teta_hub;

dev_opt_rott = dev_opt_rott + 10*(1-V_2A(end)/V_1A);
dev_opt_rotm = dev_opt_rotm + 10*(1-V_2A(end)/V_1A);
dev_opt_roth = dev_opt_roth + 10*(1-V_2A(end)/V_1A);

% INLET GEOMETRICAL ANGLE
beta_1_h_geo = beta_1_h - i_opt_roth;
beta_1_m_geo = beta_1_m - i_opt_rotm;
beta_1_t_geo = beta_1_t - i_opt_rott;

% OUTLET GEOMETRICAL ANGLE
beta_2_h_geo = beta_2_h + dev_opt_roth; 
beta_2_m_geo = beta_2_m + dev_opt_rotm;
beta_2_t_geo = beta_2_t + dev_opt_rott;

% STAGGER ANGLE
gamma_roth =  beta_2_h_geo - teta_hub/2;
gamma_rotm =  beta_2_m_geo - teta_mid/2;
gamma_rott =  beta_2_t_geo - teta_tip/2;

% ATTACK ANGLE
attack_h_design = gamma_roth - beta_1_h;
attack_m_design = gamma_rotm - beta_1_m;
attack_t_design = gamma_rott - beta_1_t;
