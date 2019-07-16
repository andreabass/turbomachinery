%% OFF-DESIGN PROBLEM (GEOMETRY)

m = 90;

% OPERATION STRATEGY: keep constant incidence angle @ rotor inlet

p1_igvinlet
p1_1D_od % Initialization (1D problem)
discretization
p1_rotin_delta
p2_rotout_od
p3_statout_od
resultsod

% ROTOR
figure(tvdrot)
subplot(3,1,3)
velt(V_1(1),W_1(1),omega*r(1),'r')
velt(V_2(1),W_2(1),omega*r(1),'r--')

subplot(3,1,2)
velt(V_1_m,W_1_m,U_m,'r')
velt(V_2_m,W_2_m,U_m,'r--')

subplot(3,1,1)
velt(V_1(end),W_1(end),omega*r(end),'r')
velt(V_2(end),W_2(end),omega*r(end),'r--')

