clear
clc

tol = 1e-6;

%% DATA
m = 100;
p_T0 = 1 .* 10.^5; % [Pa]
p_T0_t = p_T0; % [Pa]
p_T0_m = p_T0; % [Pa]
p_T0_h = p_T0; % [Pa]
T_T0 = 300; % [K]
T_T0_t = T_T0; % [K]
T_T0_m = T_T0; % [K]
T_T0_h = T_T0; % [K]
gamma = 1.4; %for air
MM = 28.84; %evaluated considering a mixture of 0.21 O2 and 0.79 N2
R_gas = 8314; % [J./kmol./K]
cp = R_gas.*gamma./(gamma-1)./MM; % [J./kgK]
R_star = R_gas ./ MM; % [J./kg./K]
beta_TT = 1.45;
V_0T_m = 0;
V_0T_t = 0; 
V_0T_h = 0;

%% ASSUMPTIONS (ISENTROPIC CASE)
D_t = 1;
lambda = 0.37;
p_0 = p_T0;
beta_IGV_m = 1.1;
V1Am_V0 = 1.11;
n = 4350;

%% CALCULATIONS
omega = n .* 2 .* pi ./ 60;
D_h = D_t.*lambda; 
b = (D_t-D_h)./2; 
D_m = (D_h+D_t)./2; 
rho_T0 = p_T0 ./ R_star ./ T_T0; 
p1_igvinlet, 
V_0 = V_0_m; 
V_1A = V_0.*V1Am_V0;
T_0 = T_0_m;
T_1_m = T_0 .* beta_IGV_m.^((1-gamma)./gamma);
V_1_m = sqrt(2.*cp.*(T_T0-T_1_m));
V_1T_m = sqrt(V_1_m.^2-V_1A.^2);
V_1T_t = V_1T_m .* D_m ./ D_t;
W_1T_t = V_1T_t - omega.*D_t./2;
W_1_t = sqrt(W_1T_t.^2 + V_1A.^2);
V_1_t = sqrt(V_1A.^2 + V_1T_t.^2);
T_1_t = T_T0 - V_1_t.^2./2./cp;
M_R1_t = W_1_t./sqrt(gamma.*R_star.*T_1_t)



