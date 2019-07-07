%%%% REFERENCE: pages 125/128 Aungier %%%

% IGV SOLIDITY
sigma_IGV_h = c_IGV / s_IGV_h;
sigma_IGV_m = c_IGV / s_IGV_m;
sigma_IGV_t = c_IGV / s_IGV_t;

% INCIDENCE ANGLE
i_IGV_h = 0;
i_IGV_m = 0;
i_IGV_t = 0;

% FLOW DEFLECTION
Dalpha_01_t = abs(alpha_0_t-alpha_1_t);
Dalpha_01_m = abs(alpha_0_m-alpha_1_m);
Dalpha_01_h = abs(alpha_0_h-alpha_1_h);

Ksh = 1;

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma_IGV_t^3/160;
p_mid = 0.914+sigma_IGV_m^3/160;
p_hub = 0.914+sigma_IGV_h^3/160;

i0_10tip = abs(alpha_0_t)^p_tip/(5+46*exp(-2.3*sigma_IGV_t))-0.1*sigma_IGV_t^3*exp((abs(alpha_0_t)-70)/4);
i0_10mid = abs(alpha_0_m)^p_mid/(5+46*exp(-2.3*sigma_IGV_m))-0.1*sigma_IGV_m^3*exp((abs(alpha_0_m)-70)/4);
i0_10hub = abs(alpha_0_h)^p_hub/(5+46*exp(-2.3*sigma_IGV_h))-0.1*sigma_IGV_h^3*exp((abs(alpha_0_h)-70)/4);

n_LIEt = 0.025*sigma_IGV_t-0.06-(abs(alpha_0_t)/90)^(1+1.2*sigma_IGV_t)/(1.5+0.43*sigma_IGV_t);
n_LIEm = 0.025*sigma_IGV_m-0.06-(abs(alpha_0_m)/90)^(1+1.2*sigma_IGV_m)/(1.5+0.43*sigma_IGV_m);
n_LIEh = 0.025*sigma_IGV_h-0.06-(abs(alpha_0_h)/90)^(1+1.2*sigma_IGV_h)/(1.5+0.43*sigma_IGV_h);

i_0_rott = Ksh*Kth*i0_10tip;
i_0_rotm = Ksh*Kth*i0_10mid;
i_0_roth = Ksh*Kth*i0_10hub;

% Deviation: dev_opt_rot = Ksh*Kth*dev0_10+m_LIE*teta

x_LIEt = abs(alpha_0_t)/100;
x_LIEm = abs(alpha_0_m)/100;
x_LIEh = abs(alpha_0_h)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;
m10_m = 0.17-0.0333*x_LIEm+0.333*x_LIEm^2;
m10_h = 0.17-0.0333*x_LIEh+0.333*x_LIEh^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;
b_m = 0.9625-0.17*x_LIEm-0.85*x_LIEm^3;
b_h = 0.9625-0.17*x_LIEh-0.85*x_LIEh^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma_IGV_t)^b_t;
m_m = m10_m/(sigma_IGV_m)^b_m;
m_h = m10_h/(sigma_IGV_h)^b_h;

dev0_10t = 0.01*sigma_IGV_t*abs(alpha_0_t)+(0.74*sigma_IGV_t^1.9+3*sigma_IGV_t)*(abs(alpha_0_t)/90)^(1.67+1.09*sigma_IGV_t);
dev0_10m = 0.01*sigma_IGV_m*abs(alpha_0_m)+(0.74*sigma_IGV_m^1.9+3*sigma_IGV_m)*(abs(alpha_0_m)/90)^(1.67+1.09*sigma_IGV_m);
dev0_10h = 0.01*sigma_IGV_h*abs(alpha_0_h)+(0.74*sigma_IGV_h^1.9+3*sigma_IGV_h)*(abs(alpha_0_h)/90)^(1.67+1.09*sigma_IGV_h);

dev_0_rott = Ksh*Kth_dev*dev0_10t;
dev_0_rotm = Ksh*Kth_dev*dev0_10m;
dev_0_roth = Ksh*Kth_dev*dev0_10h;

% CAMBER ANGLE
teta_tip_IGV = (Dalpha_01_t- i_0_rott + dev_0_rott)/(1-m_t+n_LIEt);
teta_mid_IGV = (Dalpha_01_m- i_0_rotm + dev_0_rotm)/(1-m_m+n_LIEm);
teta_hub_IGV = (Dalpha_01_h- i_0_roth + dev_0_roth)/(1-m_h+n_LIEh);

% DEVIATION ANGLE
a_c = 0.375; % Point of maximum camber NACA A4K6 (ref: Aungier pg. 70)

dev_opt_IGVh = [1 0];
while abs( (dev_opt_IGVh(end)-dev_opt_IGVh(end-1))/dev_opt_IGVh(end-1) )  > tol
    dev_opt_IGVh(end+1) = ( 0.92*a_c^2 + 0.002 * (alpha_1_h - dev_opt_IGVh(end)) ) / ( 1 - 0.002 * teta_hub_IGV / sqrt(sigma_IGV_h) ) * teta_hub_IGV / sqrt(sigma_IGV_h) + ( Ksh * Kth_dev - 1 ) * dev0_10h;
end
dev_opt_IGVh = dev_opt_IGVh(end);

dev_opt_IGVm = [1 0];
while abs( (dev_opt_IGVm(end)-dev_opt_IGVm(end-1))/dev_opt_IGVm(end-1) )  > tol
    dev_opt_IGVm(end+1) = ( 0.92*a_c^2 + 0.002 * (alpha_1_m - dev_opt_IGVm(end)) ) / ( 1 - 0.002 * teta_mid_IGV / sqrt(sigma_IGV_m) ) * teta_mid_IGV / sqrt(sigma_IGV_m) + ( Ksh * Kth_dev - 1 ) * dev0_10m;
end
dev_opt_IGVm = dev_opt_IGVm(end);

dev_opt_IGVt = [1 0];
while abs( (dev_opt_IGVt(end)-dev_opt_IGVt(end-1))/dev_opt_IGVt(end-1) )  > tol
    dev_opt_IGVt(end+1) = ( 0.92*a_c^2 + 0.002 * (alpha_1_t - dev_opt_IGVt(end)) ) / ( 1 - 0.002 * teta_tip_IGV / sqrt(sigma_IGV_t) ) * teta_tip_IGV / sqrt(sigma_IGV_t) + ( Ksh * Kth_dev - 1 ) * dev0_10t;
end
dev_opt_IGVt = dev_opt_IGVt(end);


% STAGGER ANGLE
gamma_IGVh = abs(alpha_0_h)-abs(teta_hub_IGV)/2;
gamma_IGVm = abs(alpha_0_m)-abs(teta_mid_IGV)/2;
gamma_IGVt = abs(alpha_0_t)-abs(teta_tip_IGV)/2;