clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% AXIAL COMPRESSOR DESIGN %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preliminary data and calculations

m = 100; %kg/s

%HP: Uniform temperature at IGV inlet
T_T0 = 300; %K
T_T0_t = T_T0;
T_T0_m = T_T0;
T_T0_h = T_T0;

%HP: Uniform pressure at IGV inlet
p_T0 = 1 * 10^5; %Pa
p_T0_t = p_T0;
p_T0_m = p_T0;
p_T0_h = p_T0;

beta_TT = 1.45;
gamma = 1.4; %for air
MM = 28.84; %evaluated considering a mixture of 0.21 O2 and 0.79 N2
R_gas = 8314; %J/kmol/K
cp = R_gas*gamma/(gamma-1)/MM; %J/kgK
R_star = R_gas / MM; % J/kg/K

deltaHis_TT = cp * T_T0 * beta_TT ^ ((gamma - 1)/gamma); % J/kg

%Density calculation
rho_T0 = p_T0 / R_star / T_T0; % kg/m3
rho_T0_t = p_T0_t / R_star / T_T0_t; % kg/m3
rho_T0_m = p_T0_m / R_star / T_T0_m; % kg/m3
rho_T0_h = p_T0_h / R_star / T_T0_h; % kg/m3

Q = m/rho_T0;

%% Balje diagram values
%Preliminary assumptions
D_t = 1; %m
n = 9250; %rpm

omega = n * 2 * pi / 60; %rad/s

%Calculations
D_s = D_t * (deltaHis_TT^(1/4))/sqrt(Q);
omega_s = omega * sqrt(Q)/(deltaHis_TT^(3/4));

omega_s_cordier = [10 8 7 6 5 4 3 2 1.5 1 0.9 0.8 0.7 0.6];
D_s_cordier = [0.98 1.04 1.1 1.2 1.35 1.5 1.65 1.8 2 2.2 3 3.15 3.5 4 4.75];
%
%
%
eta_TT = 0.8;

%% IGV + ROTOR INLET DESIGN

%1st subproblem

V_0A = 170;

V_0T_t = 0;
V_0T_m = 0;
V_0T_h = 0;

V_0A_t = V_0A;
V_0A_m = V_0A;
V_0A_h = V_0A;

alpha_0_t = atan(V_0T_t / V_0A_t);
alpha_0_m = atan(V_0T_m / V_0A_m);
alpha_0_h = atan(V_0T_h / V_0A_h);

V_0_t = sqrt(V_0A_t^2 + V_0T_t^2);
V_0_m = sqrt(V_0A_m^2 + V_0T_m^2);
V_0_h = sqrt(V_0A_h^2 + V_0T_h^2);

T_0_t = T_T0_t - (V_0_t^2)/(2*cp);
T_0_m = T_T0_m - (V_0_m^2)/(2*cp);
T_0_h = T_T0_h - (V_0_h^2)/(2*cp);

p_0_t = p_T0_t / (1 + (V_0_t^2)/(2 * R_star * T_0_t));
p_0_m = p_T0_m / (1 + (V_0_m^2)/(2 * R_star * T_0_m));
p_0_h = p_T0_h / (1 + (V_0_h^2)/(2 * R_star * T_0_h));

rho_0_t = p_0_t / R_star / T_0_t;
rho_0_m = p_0_m / R_star / T_0_m;
rho_0_h = p_0_h / R_star / T_0_h;

bdm = m / pi * 3 / V_0A / (rho_0_t + rho_0_h + rho_0_m);
D_h = sqrt(D_t^2 - 4*bdm);
b = (D_t - D_h) / 2;
D_m = (D_t + D_h) / 2;

%2nd subproblem
V_1A = 180;

U_1_t = omega / 2 * D_t;
U_1_m = omega / 2 * D_m;
U_1_h = omega / 2 * D_h;

phi_1_m = 0.5; %Vavra hypothesis
V_1_m = phi_1_m * U_1_m;

T_T1_t = T_T0_t;
T_T1_m = T_T0_m;
T_T1_h = T_T0_h;

V_1A_t = V_1A;
V_1A_m = V_1A;
V_1A_h = V_1A;

V_1T_m = sqrt(V_1_m^2 - V_1A_m^2);
V_1T_t = V_1T_m * D_m / D_t;
V_1T_h = V_1T_m * D_m / D_h;

alpha_1_t = atan(V_1T_t / V_1A_t);
alpha_1_m = atan(V_1T_m / V_1A_m);
alpha_1_h = atan(V_1T_h / V_1A_h);

V_1_t = sqrt(V_1A_t^2 + V_1T_t^2);
V_1_h = sqrt(V_1A_h^2 + V_1T_h^2);

W_1T_t = V_1T_t - U_1_t;
W_1T_m = V_1T_m - U_1_m;
W_1T_h = V_1T_h - U_1_h;

W_1A_t = V_1A_t;
W_1A_m = V_1A_m;
W_1A_h = V_1A_h;

W_1_t = sqrt(W_1A_t^2 + W_1T_t^2);
W_1_m = sqrt(W_1A_m^2 + W_1T_m^2);
W_1_h = sqrt(W_1A_h^2 + W_1T_h^2);

beta_1_t = atan(W_1T_t / W_1A_t);
beta_1_m = atan(W_1T_m / W_1A_m);
beta_1_h = atan(W_1T_h / W_1A_h);

%3rd subproblem

T_1_t = T_T1_t - (V_1_t^2) / (2*cp);
T_1_m = T_T1_m - (V_1_m^2) / (2*cp);
T_1_h = T_T1_h - (V_1_h^2) / (2*cp);


