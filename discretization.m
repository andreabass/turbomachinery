nodes = 1000;
rlow = linspace(D_h/2,D_m/2,nodes/2);
rhigh = linspace(D_m/2,D_t/2,nodes/2);
Dr = abs(rlow(2)-rlow(1));
c = polyfit([D_h/2 D_m/2 D_t/2],[alpha_1_h_geo alpha_1_m_geo alpha_1_t_geo],length([D_h/2 D_m/2 D_t/2])-1);
alpha_1_geo_low = polyval(c,rlow);
alpha_1_geo_high = polyval(c,rhigh);
