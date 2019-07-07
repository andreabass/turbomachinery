
Dalpha_23_t = abs(alpha_2_t-alpha_3_t);
Dalpha_23_m = abs(alpha_2_m-alpha_3_m);
Dalpha_23_h = abs(alpha_2_h-alpha_3_h);

% Incidence: i_opt_stat = Ksh*Kth*i0_10+n_LIE*teta
Ksh = 1; % for NACA-65

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma_S_t^3/160;
p_mid = 0.914+sigma_S_m^3/160;
p_hub = 0.914+sigma_S_h^3/160;

i0_10tip = abs(alpha_2_t)^p_tip/(5+46*exp(-2.3*sigma_S_t))-0.1*sigma_S_t^3*exp((abs(alpha_2_t)-70)/4);
i0_10mid = abs(alpha_2_m)^p_mid/(5+46*exp(-2.3*sigma_S_m))-0.1*sigma_S_m^3*exp((abs(alpha_2_m)-70)/4);
i0_10hub = abs(alpha_2_h)^p_hub/(5+46*exp(-2.3*sigma_S_h))-0.1*sigma_S_h^3*exp((abs(alpha_2_h)-70)/4);

n_LIEt = 0.025*sigma_S_t-0.06-(abs(alpha_2_t)/90)^(1+1.2*sigma_S_t)/(1.5+0.43*sigma_S_t);
n_LIEm = 0.025*sigma_S_m-0.06-(abs(alpha_2_m)/90)^(1+1.2*sigma_S_m)/(1.5+0.43*sigma_S_m);
n_LIEh = 0.025*sigma_S_h-0.06-(abs(alpha_2_h)/90)^(1+1.2*sigma_S_h)/(1.5+0.43*sigma_S_h);

i_0_statt = Ksh*Kth*i0_10tip;
i_0_statm = Ksh*Kth*i0_10mid;
i_0_stath = Ksh*Kth*i0_10hub;

% Deviation: dev_opt_stat = Ksh*Kth*dev0_10+m_LIE*teta

x_LIEt = abs(alpha_2_t)/100;
x_LIEm = abs(alpha_2_m)/100;
x_LIEh = abs(alpha_2_h)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;
m10_m = 0.17-0.0333*x_LIEm+0.333*x_LIEm^2;
m10_h = 0.17-0.0333*x_LIEh+0.333*x_LIEh^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;
b_m = 0.9625-0.17*x_LIEm-0.85*x_LIEm^3;
b_h = 0.9625-0.17*x_LIEh-0.85*x_LIEh^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma_S_t)^b_t;
m_m = m10_m/(sigma_S_m)^b_m;
m_h = m10_h/(sigma_S_h)^b_h;

dev0_10t = 0.01*sigma_S_t*abs(alpha_2_t)+(0.74*sigma_S_t^1.9+3*sigma_S_t)*(abs(alpha_2_t)/90)^(1.67+1.09*sigma_S_t);
dev0_10m = 0.01*sigma_S_m*abs(alpha_2_m)+(0.74*sigma_S_m^1.9+3*sigma_S_m)*(abs(alpha_2_m)/90)^(1.67+1.09*sigma_S_m);
dev0_10h = 0.01*sigma_S_h*abs(alpha_2_h)+(0.74*sigma_S_h^1.9+3*sigma_S_h)*(abs(alpha_2_h)/90)^(1.67+1.09*sigma_S_h);

dev_0_statt = Ksh*Kth_dev*dev0_10t;
dev_0_statm = Ksh*Kth_dev*dev0_10m;
dev_0_stath = Ksh*Kth_dev*dev0_10h;

% CAMBER ANGLE
teta_tip_stat = (Dalpha_23_t- i_0_statt + dev_0_statt)/(1-m_t+n_LIEt);
teta_mid_stat = (Dalpha_23_m- i_0_statm + dev_0_statm)/(1-m_m+n_LIEm);
teta_hub_stat = (Dalpha_23_h- i_0_stath + dev_0_stath)/(1-m_h+n_LIEh);

% INCIDENCE ANGLE
i_opt_statt = Ksh*Kth*i0_10tip+n_LIEt*teta_tip_stat;
i_opt_statm = Ksh*Kth*i0_10mid+n_LIEm*teta_mid_stat;
i_opt_stath = Ksh*Kth*i0_10hub+n_LIEh*teta_hub_stat;

% DEVIATION ANGLE
dev_opt_statt = Ksh*Kth_dev*dev0_10t+m_t*teta_tip_stat;
dev_opt_statm = Ksh*Kth_dev*dev0_10m+m_m*teta_mid_stat;
dev_opt_stath = Ksh*Kth_dev*dev0_10h+m_h*teta_hub_stat;

% STAGGER ANGLE
gamma_stath = abs(alpha_2_h)-abs(teta_hub_stat)/2+abs(i_opt_stath);
gamma_statm = abs(alpha_2_m)-abs(teta_mid_stat)/2+abs(i_opt_statm);
gamma_statt = abs(alpha_2_t)-abs(teta_tip_stat)/2+abs(i_opt_statt);
