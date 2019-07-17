function i_min = i_min(beta_in,beta_out,th_c,sigma)

Dbeta = abs(beta_in-beta_out);

Ksh_i = 0.7; % for NACA-65

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma^3/160;

i0_10tip = abs(beta_in)^p_tip/(5+46*exp(-2.3*sigma))-0.1*sigma^3*exp((abs(beta_in)-70)/4);

n_LIEt = 0.025*sigma-0.06-(abs(beta_in)/90)^(1+1.2*sigma)/(1.5+0.43*sigma);


i_0_rott = Ksh_i*Kth*i0_10tip;


% Deviation: dev_opt_rot = Ksh*Kth*dev0_10+m_LIE*teta

Ksh_d = 0.75;

x_LIEt = abs(beta_in)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma)^b_t;

dev0_10t = 0.01*sigma*abs(beta_in)+(0.74*sigma^1.9+3*sigma)*(abs(beta_in)/90)^(1.67+1.09*sigma);

dev_0_rott = Ksh_d*Kth_dev*dev0_10t;

% CAMBER ANGLE
teta_tip = (Dbeta- i_0_rott + dev_0_rott)/(1-m_t+n_LIEt);

% INCIDENCE ANGLE
i_min = Ksh_i*Kth*i0_10tip+n_LIEt*teta_tip;

end

