nodes = 100;
rlow  = linspace(D_m/2,D_h/2,nodes/2);
rhigh = linspace(D_m/2,D_t/2,nodes/2);
r     = [rlow(end:-1:1) rhigh(2:end) ];
Drhigh = rhigh(2)-rhigh(1);
Drlow = rlow(2)-rlow(1);
Dr = Drhigh;

% alpha_1_geo
c = polyfit([D_h/2 D_m/2 D_t/2],[alpha_1_h_geo alpha_1_m_geo alpha_1_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
alpha_1_geo_low = polyval(c,rlow);
alpha_1_geo_high = polyval(c,rhigh);

% dev_opt_rot
c = polyfit([D_h/2 D_m/2 D_t/2],[dev_opt_roth dev_opt_rotm dev_opt_rott],length([D_h/2 D_m/2 D_t/2])-1);
dev_opt_rot_low = polyval(c,rlow);
dev_opt_rot_high = polyval(c,rhigh);

% i_opt_rot
c = polyfit([D_h/2 D_m/2 D_t/2],[i_opt_roth i_opt_rotm i_opt_rott],length([D_h/2 D_m/2 D_t/2])-1);
i_opt_rot_low = polyval(c,rlow);
i_opt_rot_high = polyval(c,rhigh);

% beta_1_geo
c = polyfit([D_h/2 D_m/2 D_t/2],[beta_1_h_geo beta_1_m_geo beta_1_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
beta_1_geo_low = polyval(c,rlow);
beta_1_geo_high = polyval(c,rhigh);

% beta_2_geo
c = polyfit([D_h/2 D_m/2 D_t/2],[beta_2_h_geo beta_2_m_geo beta_2_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
beta_2_geo_low = polyval(c,rlow);
beta_2_geo_high = polyval(c,rhigh);

% gamma_rot
c = polyfit([D_h/2 D_m/2 D_t/2],[gamma_roth gamma_rotm gamma_rott],length([D_h/2 D_m/2 D_t/2])-1);
gamma_rot_low = polyval(c,rlow);
gamma_rot_high = polyval(c,rhigh);

% attack_design (rotor)
c = polyfit([D_h/2 D_m/2 D_t/2],[attack_h_design attack_m_design attack_t_design],length([D_h/2 D_m/2 D_t/2])-1);
attack_design_low = polyval(c,rlow);
attack_design_high = polyval(c,rhigh);

% teta (rotor)
c = polyfit([D_h/2 D_m/2 D_t/2],[teta_hub teta_mid teta_tip],length([D_h/2 D_m/2 D_t/2])-1);
teta_low = polyval(c,rlow);
teta_high = polyval(c,rhigh);

%alpha_2_geo
c = polyfit([D_h/2 D_m/2 D_t/2],[alpha_2_h_geo alpha_2_m_geo alpha_2_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
alpha_2_geo_low = polyval(c,rlow);
alpha_2_geo_high = polyval(c,rhigh);

%alpha_3_geo
c = polyfit([D_h/2 D_m/2 D_t/2],[alpha_3_h_geo alpha_3_m_geo alpha_3_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
alpha_3_geo_low = polyval(c,rlow);
alpha_3_geo_high = polyval(c,rhigh);

% dev_opt_stat
c = polyfit([D_h/2 D_m/2 D_t/2],[dev_opt_stath dev_opt_statm dev_opt_statt],length([D_h/2 D_m/2 D_t/2])-1);
dev_opt_stat_low = polyval(c,rlow);
dev_opt_stat_high = polyval(c,rhigh);

% i_opt_stat
c = polyfit([D_h/2 D_m/2 D_t/2],[i_opt_stath i_opt_statm i_opt_statt],length([D_h/2 D_m/2 D_t/2])-1);
i_opt_stat_low = polyval(c,rlow);
i_opt_stat_high = polyval(c,rhigh);

% gamma_stat
c = polyfit([D_h/2 D_m/2 D_t/2],[gamma_stath gamma_statm gamma_statt],length([D_h/2 D_m/2 D_t/2])-1);
gamma_stat_low = polyval(c,rlow);
gamma_stat_high = polyval(c,rhigh);

% attack_design (stator)
c = polyfit([D_h/2 D_m/2 D_t/2],[attack_S_h_design attack_S_m_design attack_S_t_design],length([D_h/2 D_m/2 D_t/2])-1);
attack_S_design_low = polyval(c,rlow);
attack_S_design_high = polyval(c,rhigh);

% teta (stator)
c = polyfit([D_h/2 D_m/2 D_t/2],[teta_hub_stat teta_mid_stat teta_tip_stat],length([D_h/2 D_m/2 D_t/2])-1);
teta_S_low = polyval(c,rlow);
teta_S_high = polyval(c,rhigh);