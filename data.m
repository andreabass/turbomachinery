%% Preliminary data and calculations

m = 100; % [kg/s]

%Data: uniform total temperature at IGV inlet
T_T0_m = 300;    % [K] (eq. 1)
T_T0   = T_T0_m; % (definition)
T_T0_t = T_T0;   % (eq.2)
T_T0_h = T_T0;   % (eq.3)

%Data: uniform total pressure at IGV inlet
p_T0_m = 1 * 10^5; % [Pa] (eq. 4)
p_T0   = p_T0_m;   % (definition)
p_T0_t = p_T0;     % (eq. 5)
p_T0_h = p_T0;     % (eq. 6)

beta_TT = 1.45;
gamma = 1.4; %for air
MM = 28.84; %evaluated considering a mixture of 0.21 O2 and 0.79 N2
R_gas = 8314; % [J/kmol/K]
cp = R_gas*gamma/(gamma-1)/MM; % [J/kgK]
R_star = R_gas / MM; % [J/kg/K]
mu = 1.81*10^-5; % [Pas]

deltaHis_TT = cp * T_T0 * (beta_TT ^ ((gamma - 1)/gamma)-1); % [J/kg]

%Density calculation
rho_T0_m = p_T0_m / R_star / T_T0_m; % [kg/m3] (eq. 7)
rho_T0_t = p_T0_t / R_star / T_T0_t; % [kg/m3] (eq. 8)
rho_T0_h = p_T0_h / R_star / T_T0_h; % [kg/m3] (eq. 9)
rho_T0   = rho_T0_m;                 % (definition)