%------------------------------------------------------------------------
%------------------------------------------------------------------------
%                   AXIAL COMPRESSOR DESIGN
%------------------------------------------------------------------------
%------------------------------------------------------------------------

clc
clear

%~~~~~ Assumption ~~~~~
% - Free vortex design for IGV, rotor and stator
% - Axial inlet alpha_0 = 0∞

% DATA

m = 100; %kg/s
T_t0 = 300; %K
P_t0 = 1; %bar
beta_tt = 1.45;
gamma = 1.4; %assumed; mixture of diatomic species
MMa = 28.84; %evaluated considering 0.21 O2 and 0.79 N2
Rgas = 8314; %J/kmol/K
cpa = Rgas*gamma/(gamma-1)/MMa; %J/kgK

%------------------------------------------------------------------------

% BALJE APPROACH

Dhs_tt = cpa*T_t0*(beta_tt^((gamma-1)/gamma)-1); %J/kg
rho_t0 = P_t0*10^5/(Rgas/MMa*T_t0); %kg/m3

QBalje = m/rho_t0;

%~~~~~ Assumption ~~~~~
Dtip = 1; % limited at the maximum of 1 m
n = 9250; % rpm
omega = 2*pi*n/60; 

% At first we have chosen a value for omega_s = 5 (the value falling
% approximately on the Chordier line for a fixed D_s =
% Dtip*Dhs_tt^(1/4)/sqrt(QBalje); in this case peripheral velocity in this
% case turns to be too high (~640 m/s). In order to reduce this value, we
% have accepted a slightly reduction (see Balje chart) in efficiency,
% fixing the numer of round per minute to 6000 rpm (gear box needed);
% eventually the corresponding value for omega_s has been calculated
%NOTE: if n too low, some problems occurs in 'while' cycle for b_1, since
%decreasing n-->omega decrease-->b_1vavra increase: it can increase so much
%that Dhub becomes negative (not possible in reality).
%~~~~~~~~~~~~~~~~~~~~~~~

D_s = Dtip*Dhs_tt^(1/4)/sqrt(QBalje);
omega_s = omega*sqrt(QBalje)/Dhs_tt^(3/4);

eta_tt = [0 0.8]; % from Balje ----> must be checked at the end
lambda = 0.5; % from Balje

% Overall iterative process: check on eta_TT (initially
% assumed from Balje) and P_3av, that must be equal to the required by the
% costumer (P_t3av = beta_tt*P_t0).

P_t3_av = 2; % first guess to enter the cycle
P_t3m = 2;
tol = 1e-4;
tol0 = 1e-3;
g = 1;

% VAVRA

phi = 0.5; % flow coefficient (v_ax/u_mid) at midspan
X = 0.5; % reaction degree at midspan

 while abs(eta_tt(g+1)-eta_tt(g))>tol || abs(P_t3m(end)-beta_tt*P_t0)>tol0

Leul = Dhs_tt/eta_tt(g+1); % J/kg

rho_1m = [0 rho_t0]; % value to cbe changed in an iterative process from IGV
j = 1;
tol1 = 1e-4;

b_1 = 0; % value set to enter the 'while' cycle
b_1vavr = 1; % value set to enter the 'while' cycle

while abs(rho_1m(j+1)-rho_1m(j))>tol1

Dhub = lambda*Dtip; 
i = 1;
tol2 = 1e-4;

while abs(b_1vavr-b_1)>tol2 % iterative cycle to get same value for b1Vavra and b1
b_1 = (Dtip-Dhub)/2;
Dmid =(Dtip+Dhub)/2;
u_mid = omega*Dmid/2;
v_1ax = phi*u_mid; % constant along the span (free vortex method)
b_1vavr = m/(rho_1m(j+1)*pi*Dmid*v_1ax);
Dhub = Dtip-2*b_1vavr;
lambda = Dhub/Dtip;
end

% x0 = [100 100]; % intial values for v1t and v2t
% Sys = sistema
% [v_1tm v_2tm] = fsolve(@sys, x0)

syms v1t v2t
sol1 = solve((v1t+v2t)/(2*u_mid)+X-1==0, v2t-v1t-Leul/u_mid==0,[v1t,v2t]);
v_1tm = double(sol1.v1t); % to convert the solution from symbolic to numeric
v_2tm = double(sol1.v2t); % to convert the solution from symbolic to numeric


%------------------------------------------------------------------------
%                                IGV
%------------------------------------------------------------------------

v_1m = sqrt(v_1ax^2+v_1tm^2);
alpha_1m = atand(v_1tm/v_1ax);
T_t1 = T_t0; % IGV is a statoric blade (no work exchanged)
T_1m = T_t1-v_1m^2/(2*cpa); % K
M_1m = v_1m/sqrt(gamma*Rgas/MMa*(T_1m));

%~~~~~ Assumption ~~~~~
alpha_0 = 0;
Xmid_IGV = 0.5;
chord = 0.04; %m
%~~~~~~~~~~~~~~~~~~~~~~

% From AINLEY-MATHIESON correlation
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5
% X_Anl,Math = s/c-s/c_min 0 since for us s/c=s/c_min
% B,C,n not considered for the same reason as above

alpha_2primem = 90-alpha_1m;

if alpha_2primem <= 30
    s_c_min = 0.46+alpha_2primem/77; % ratio s/c which minimize the profile losses
else 
    s_c_min = 0.614+alpha_2primem/130; % ratio s/c which minimize the profile losses
end 

if alpha_2primem <= 27
    A = 0.025+(27-alpha_2primem)/530;
else
    A = 0.025+(27-alpha_2primem)/3850;
end

Yp1 = A;

mua = 1.81*10^-5; %Pas 
Re_ref = 2*10^5;
Re_m = rho_1m(j+1)*v_1m*chord/mua; % Re at midspan
alpha_av01 = atand((tand(alpha_0)+tand(alpha_1m))/2);
cL = 2*s_c_min*(abs(tand(alpha_1m)-tand(alpha_0)))*cosd(alpha_av01);

Yp1_Re = Yp1*(Re_ref/Re_m)^0.2; % Re correction for Yp1
Ysec = chord/b_1vavr*(0.0334*cosd(alpha_1m)/cosd(alpha_0))*(cL/s_c_min)^2*(cosd(alpha_1m))^2/(cosd(alpha_av01))^3; % since it's a statoric blade no clearance losses
Ytot = Yp1_Re+Ysec;

syms Pt1 P1
sol2 = solve(Ytot-(P_t0-Pt1)/(Pt1-P1)==0, Pt1/P1-(1+(gamma-1)/2*M_1m^2)^(gamma/(gamma-1))==0,[Pt1,P1]);
P_t1m = double(sol2.Pt1); % to convert the solution from symbolic to numeric
P_1m = double(sol2.P1); % to convert the solution from symbolic to numeric

rho_1m = [rho_1m P_1m*10^5/(Rgas/MMa*T_1m)];

b_1vavr = m/(rho_1m(end)*pi*Dmid*v_1ax);
Dhub = Dtip-2*b_1vavr;
b_1 = (Dtip-Dhub)/2;

j=j+1;

end

% Inlet IGV: iterative cycle to evaluate inlet IGV condition

b_0 = b_1vavr; % constant height of the blade
rho_0 = [0 rho_t0]; % first guess
z = 1;
tol3 = 1e-4;

while abs(rho_0(z+1)-rho_0(z))>tol3
v_0 = m/(rho_0(z+1)*pi*Dmid*b_0); % constant along the span (axial)
T_0 = T_t0-v_0^2/(2*cpa);
M_0 = v_0/sqrt(gamma*Rgas/MMa*T_0);
P_0 = P_t0/((1+(gamma-1)/2*M_0^2)^(gamma/(gamma-1)));

rho_0 = [rho_0 P_0*10^5/(Rgas/MMa*T_0)];

z = z+1;
end

% number of IGV blades
s_IGVm = chord*s_c_min; % pitch at midspan
Nblade_IGV = ceil(pi*Dmid/s_IGVm); % neglecting the thickness

%------------------------------------------------------------------------
%                            ROTOR INLET
%------------------------------------------------------------------------

%~~~~~ Assumption ~~~~~
v_2ax = v_1ax;
w_1ax = v_1ax;
w_2ax = v_2ax;
% all evaluations here proposed are performed at midspan
%~~~~~~~~~~~~~~~~~~~~~~

% Check on HOWELL correlation: check performed in a graphical way (more
% accurate correlation will be used further on)
v_2m = sqrt(v_2ax^2+v_2tm^2);
alpha_2m = atand(v_2tm/v_2ax);
Dalpha_m12 = abs(alpha_2m-alpha_1m);
w_1tm = v_1tm-u_mid;
w_1m = sqrt(w_1tm^2+w_1ax^2);
beta_1m = atand(w_1tm/w_1ax);
w_2tm = v_2tm-u_mid;
beta_2m = atand(w_2tm/w_2ax);

Dbeta_m12 = abs(beta_2m-beta_1m);

Mr_1m = w_1m/sqrt(gamma*Rgas/MMa*T_1m);

% at first iteration starting with n = 8500 rpm we have obtained results in
% term of HOWELL not that satisfynig, obtaining a condition too much
% loaded. Keeping constant Leul increasing the peripheral velocity we can
% obtain a reduction of Dvt and hence of the blade loading. The value
% chosen for n = 9000. For such value the point on HOWELL chart falls right
% on the line of optimal deflection (considering psi = phi =1).

%-------------------------------------------------------------------------

% VELOCITY TRAINGLES along the span according to free vortex criterion (
% v_t*R = const)

% section (1): HUB
u_hub = omega*Dhub/2;
v_1th = v_1tm*Dmid/Dhub;
v_1h = sqrt(v_1th^2+v_1ax^2);
w_1th = v_1th-u_hub;
w_1h = sqrt(w_1th^2+w_1ax^2);

alpha_1h = atand(v_1th/v_1ax);
beta_1h = atand(w_1th/w_1ax);

% thermodynamic quantities (losses according to AINLEY_MATHIESON)

T_1h = T_t1-v_1h^2/(2*cpa);
M_1h = v_1h/sqrt(gamma*Rgas/MMa*T_1h);
Mr_1h = w_1h/sqrt(gamma*Rgas/MMa*T_1h);

% From AINLEY-MATHIESON correlation (applied to IGV)
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5

alpha_2primeh = 90-alpha_1h;

C = 0.08*((alpha_2primeh/30)^2-1);
n_AM = 1+alpha_2primeh/30;
s_c_h = pi*Dhub/(Nblade_IGV*chord);

if alpha_2primeh <= 27
    A = 0.025+(27-alpha_2primeh)/530;
else
    A = 0.025+(27-alpha_2primeh)/3850;
end

if alpha_2primeh <= 30
    s_c_minh = 0.46+alpha_2primeh/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2primeh/1640;
    X_AM = s_c_h-s_c_minh;
    
    Yp1h = A+B*X_AM^2+C*X_AM^3; % proffile losses correlation
else 
    s_c_minh = 0.614+alpha_2primeh/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = s_c_h-s_c_minh;
    
    Yp1h = A+B*(abs(X_AM))^n_AM;
end 

% iterative process for rho1 at hub
rho_1h = [0 rho_1m(end)]; % as first guess the value of rho1 at midspan has been considered
tol4 = 1e-4;
k = 1;

while abs(rho_1h(k+1)-rho_1h(k))>tol4
Re_h = rho_1h(k+1)*v_1h*chord/mua; % Re at midspan
alpha_av01 = atand((tand(alpha_0)+tand(alpha_1h))/2);
cL = 2*s_c_h*(abs(tand(alpha_1h)-tand(alpha_0)))*cosd(alpha_av01);

Yp1h_Re = Yp1h*(Re_ref/Re_h)^0.2; % Re correction for Yp1
Ysech = chord/b_1vavr*(0.0334*cosd(alpha_1h)/cosd(alpha_0))*(cL/s_c_h)^2*(cosd(alpha_1h))^2/(cosd(alpha_av01))^3; % since it's a statoric blade no clearance losses
Ytoth = Yp1h_Re+Ysech;

syms Pt1h P1h
sol3 = solve(Ytoth-(P_t0-Pt1h)/(Pt1h-P1h)==0, Pt1h/P1h-(1+(gamma-1)/2*M_1h^2)^(gamma/(gamma-1))==0,[Pt1h,P1h]);
P_t1h = double(sol3.Pt1h); % to convert the solution from symbolic to numeric
P_1h = double(sol3.P1h); % to convert the solution from symbolic to numeric

rho_1h = [rho_1h P_1h*10^5/(Rgas/MMa*T_1h)];

k = k+1;
end

%-------------------------------------------------------------------------
% section (1): TIP
u_tip = omega*Dtip/2;
v_1tt = v_1tm*Dmid/Dtip;
v_1t = sqrt(v_1tt^2+v_1ax^2);
w_1tt = v_1tt-u_tip;
w_1t = sqrt(w_1tt^2+w_1ax^2);

alpha_1t = atand(v_1tt/v_1ax);
beta_1t = atand(w_1tt/w_1ax);

% thermodynamic quantities (losses according to AINLEY_MATHIESON)

T_1t = T_t1-v_1t^2/(2*cpa);
M_1t = v_1t/sqrt(gamma*Rgas/MMa*T_1t);
Mr_1t = w_1t/sqrt(gamma*Rgas/MMa*T_1t);

% From AINLEY-MATHIESON correlation (applied to IGV)
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5

alpha_2primet = 90-alpha_1h;

C = 0.08*((alpha_2primet/30)^2-1);
n_AM = 1+alpha_2primet/30;
s_c_t = pi*Dtip/(Nblade_IGV*chord);

if alpha_2primet <= 27
    A = 0.025+(27-alpha_2primet)/530;
else
    A = 0.025+(27-alpha_2primet)/3850;
end

if alpha_2primet <= 30
    s_c_mint = 0.46+alpha_2primet/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2primet/1640;
    X_AM = s_c_t-s_c_mint;
    
    Yp1t = A+B*X_AM^2+C*X_AM^3; % proffile losses correlation
else 
    s_c_mint = 0.614+alpha_2primet/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = s_c_t-s_c_mint;
    
    Yp1t = A+B*(abs(X_AM))^n_AM;
end 

% iterative process for rho1 at hub
rho_1t = [0 rho_1m(end)]; % as first guess the value of rho1 at midspan has been considered
tol4 = 1e-4;
k = 1;

while abs(rho_1t(k+1)-rho_1t(k))>tol4
Re_t = rho_1t(k+1)*v_1t*chord/mua; % Re at midspan
alpha_av01 = atand((tand(alpha_0)+tand(alpha_1t))/2);
cL = 2*s_c_t*(abs(tand(alpha_1t)-tand(alpha_0)))*cosd(alpha_av01);

Yp1t_Re = Yp1t*(Re_ref/Re_t)^0.2; % Re correction for Yp1
Ysect = chord/b_1vavr*(0.0334*cosd(alpha_1t)/cosd(alpha_0))*(cL/s_c_t)^2*(cosd(alpha_1t))^2/(cosd(alpha_av01))^3; % since it's a statoric blade no clearance losses
Ytott = Yp1t_Re+Ysect;

syms Pt1t P1t
sol4 = solve(Ytott-(P_t0-Pt1t)/(Pt1t-P1t)==0, Pt1t/P1t-(1+(gamma-1)/2*M_1t^2)^(gamma/(gamma-1))==0,[Pt1t,P1t]);
P_t1t = double(sol4.Pt1t); % to convert the solution from symbolic to numeric
P_1t = double(sol4.P1t); % to convert the solution from symbolic to numeric

rho_1t = [rho_1t P_1t*10^5/(Rgas/MMa*T_1t)];

k = k+1;
end

% IGV efficiency (eta_IGV_SS)

T_1tis = T_0*(P_1t/P_0)^((gamma-1)/gamma);
T_1mis = T_0*(P_1m/P_0)^((gamma-1)/gamma);
T_1his = T_0*(P_1h/P_0)^((gamma-1)/gamma);

eta_IGV_tip = (T_0-T_1t)/(T_0-T_1tis);
eta_IGV_mid = (T_0-T_1m)/(T_0-T_1mis);
eta_IGV_hub = (T_0-T_1h)/(T_0-T_1his);

%------------------------------------------------------------------------
%                            ROTOR OUTLET
%------------------------------------------------------------------------

b_2 = b_1; % assumption: constant blade height between inlet and outlet

v_2ax=[0 v_1ax];
o = 1;
tol5 = 1e-4;

eta_rot_tip = [0 eta_tt(end)];
eta_rot_mid = [0 eta_tt(end)];
eta_rot_hub = [0 eta_tt(end)];
f = 1;

% Iterative cycle on the rotor in order to find  the real rotor efficiency
% (convergence also on p_2 value at hub, mid, tip

while abs(eta_rot_tip(f+1)-eta_rot_tip(f))>tol5 || abs(eta_rot_mid(f+1)-eta_rot_mid(f))>tol5 || abs(eta_rot_hub(f+1)-eta_rot_hub(f))>tol5

% Iterative cycle considering that up to now we have considered vax1=vax2
% but this is not possible since in axial machine vax is not constant 
% v_2ax=[0 v1_ax]
  
while abs(v_2ax(o+1)-v_2ax(o))>tol5
    
% section (2): HUB

v_2th = v_2tm*Dmid/Dhub;
v_2h = sqrt(v_2th^2+v_2ax(o+1)^2);
w_2th = v_2th-u_hub;
w_2h = sqrt(w_2th^2+w_2ax^2);

alpha_2h = atand(v_2th/v_2ax(o+1));
beta_2h = atand(w_2th/w_2ax);

X_rh = (w_1ax^2+w_1th^2-w_2ax^2-w_2th^2)/(2*Leul); % remember to define w2ax in each cycle
T_2h = T_1h+X_rh*Leul/cpa;
T_2his = eta_rot_hub(f+1)*(T_2h-T_1h)+T_1h; % assuming eta_rotor = eta_tt since degree of reaction is around 0.5 (quite strong hypotesis since at the hub/tip degree of reaction is different from 0.5)
P_2h = P_1h*(T_2his/T_1h)^(gamma/(gamma-1));
rho_2h = P_2h*10^5/(Rgas/MMa*T_2h);

% section (2): MID

v_2m = sqrt(v_2ax(o+1)^2+v_2tm^2);
w_2m = sqrt(w_2ax^2+w_2tm^2);

alpha_2m = atand(v_2tm/v_2ax(o+1));
beta_2m = atand(w_2tm/w_2ax);

X_rm = (w_1ax^2+w_1tm^2-w_2ax^2-w_2tm^2)/(2*Leul); % remember to define w2ax in each cycle
T_2m = T_1m+X_rm*Leul/cpa;
T_2mis = eta_rot_mid(f+1)*(T_2m-T_1m)+T_1m; % assuming eta_rotor = eta_tt since degree of reaction is around 0.5
P_2m = P_1m*(T_2mis/T_1m)^(gamma/(gamma-1));
rho_2m = P_2m*10^5/(Rgas/MMa*T_2m);

% section (2): TIP

v_2tt = v_2tm*Dmid/Dtip;
v_2t = sqrt(v_2tt^2+v_2ax(o+1)^2);
w_2tt = v_2tt-u_tip;
w_2t = sqrt(w_2tt^2+w_2ax^2);

alpha_2t = atand(v_2tt/v_2ax(o+1));
beta_2t = atand(w_2tt/w_2ax);

X_rt = (w_1ax^2+w_1tt^2-w_2ax^2-w_2tt^2)/(2*Leul); % remember to define w2ax in each cycle
T_2t = T_1t+X_rt*Leul/cpa;
T_2tis = eta_rot_tip(f+1)*(T_2t-T_1t)+T_1t; % assuming eta_rotor = eta_tt since degree of reaction is around 0.5 (quite strong hypotesis since at the hub/tip degree of reaction is different from 0.5)
P_2t = P_1t*(T_2tis/T_1t)^(gamma/(gamma-1));
rho_2t = P_2t*10^5/(Rgas/MMa*T_2t);


rho_av_hmt = (rho_2h+rho_2m+rho_2t)/3;
v_2ax = [v_2ax m/(rho_av_hmt*pi*Dmid*b_2)];

w_2ax = v_2ax(end);

o = o+1;
end

% Check on critcal Mach number
M_2h = v_2h/sqrt(gamma*Rgas/MMa*T_2h);
M_2m = v_2m/sqrt(gamma*Rgas/MMa*T_2m);
M_2t = v_2t/sqrt(gamma*Rgas/MMa*T_2t);
Mr_2h = w_2h/sqrt(gamma*Rgas/MMa*T_2h);
Mr_2m = w_2m/sqrt(gamma*Rgas/MMa*T_2m);
Mr_2t = w_2t/sqrt(gamma*Rgas/MMa*T_2t);

% Check on HOWELL correlation: check performed in a graphical way (more
% accurate correlation will be used further on)
Dbeta_h12 = abs(beta_2h-beta_1h);
Dbeta_m12 = abs(beta_2m-beta_1m);
Dbeta_t12 = abs(beta_2t-beta_1t);

% Evaluation of the minimum chord which allows to have a Re = Re_refHOWE =
% 3*10^5
Re_refHOWE = 3*10^5;


w_av12t = (w_1t+w_2t)/2;
w_av12m = (w_1m+w_2m)/2;
w_av12h = (w_1h+w_2h)/2;

rho_av12t = (rho_1t(end)+rho_2t)/2;
rho_av12m = (rho_1m(end)+rho_2m)/2;
rho_av12h = (rho_1h(end)+rho_2h)/2;

chord_mint = Re_refHOWE*mua/(rho_av12t*w_av12t);
chord_minm = Re_refHOWE*mua/(rho_av12m*w_av12m);
chord_minh = Re_refHOWE*mua/(rho_av12h*w_av12h);

chord_rot = 0.1; % assumed

sigma_rm = 0.9; % (at mid) assumed considering the results previously obtained with HOWELL we have noticed that the blade at midspan should be more loaded hence sigma should decrease
s_rotm = chord_rot/sigma_rm;

Nblade_rot = ceil(pi*Dmid/s_rotm);

sigma_rt = chord_rot/(pi*Dtip/Nblade_rot);
sigma_rh = chord_rot/(pi*Dhub/Nblade_rot);

% LIEBLEIN 

th_c = 0.08; % percentage thickness (assumed)

% Incidence: i_opt_rot = Ksh*Kth*i0_10+n_LIE*teta
Ksh = 1; % for NACA-65

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma_rt^3/160;
p_mid = 0.914+sigma_rm^3/160;
p_hub = 0.914+sigma_rh^3/160;

i0_10tip = abs(beta_1t)^p_tip/(5+46*exp(-2.3*sigma_rt))-0.1*sigma_rt^3*exp((abs(beta_1t)-70)/4);
i0_10mid = abs(beta_1m)^p_mid/(5+46*exp(-2.3*sigma_rm))-0.1*sigma_rm^3*exp((abs(beta_1m)-70)/4);
i0_10hub = abs(beta_1h)^p_hub/(5+46*exp(-2.3*sigma_rh))-0.1*sigma_rh^3*exp((abs(beta_1h)-70)/4);

n_LIEt = 0.025*sigma_rt-0.06-(abs(beta_1t)/90)^(1+1.2*sigma_rt)/(1.5+0.43*sigma_rt);
n_LIEm = 0.025*sigma_rm-0.06-(abs(beta_1m)/90)^(1+1.2*sigma_rm)/(1.5+0.43*sigma_rm);
n_LIEh = 0.025*sigma_rh-0.06-(abs(beta_1h)/90)^(1+1.2*sigma_rh)/(1.5+0.43*sigma_rh);

i_0_rott = Ksh*Kth*i0_10tip;
i_0_rotm = Ksh*Kth*i0_10mid;
i_0_roth = Ksh*Kth*i0_10hub;

% Deviation: dev_opt_rot = Ksh*Kth*dev0_10+m_LIE*teta

x_LIEt = abs(beta_1t)/100;
x_LIEm = abs(beta_1m)/100;
x_LIEh = abs(beta_1h)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;
m10_m = 0.17-0.0333*x_LIEm+0.333*x_LIEm^2;
m10_h = 0.17-0.0333*x_LIEh+0.333*x_LIEh^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;
b_m = 0.9625-0.17*x_LIEm-0.85*x_LIEm^3;
b_h = 0.9625-0.17*x_LIEh-0.85*x_LIEh^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma_rt)^b_t;
m_m = m10_m/(sigma_rm)^b_m;
m_h = m10_h/(sigma_rh)^b_h;

dev0_10t = 0.01*sigma_rt*abs(beta_1t)+(0.74*sigma_rt^1.9+3*sigma_rt)*(abs(beta_1t)/90)^(1.67+1.09*sigma_rt);
dev0_10m = 0.01*sigma_rm*abs(beta_1m)+(0.74*sigma_rm^1.9+3*sigma_rm)*(abs(beta_1m)/90)^(1.67+1.09*sigma_rm);
dev0_10h = 0.01*sigma_rh*abs(beta_1h)+(0.74*sigma_rh^1.9+3*sigma_rh)*(abs(beta_1h)/90)^(1.67+1.09*sigma_rh);

dev_0_rott = Ksh*Kth_dev*dev0_10t;
dev_0_rotm = Ksh*Kth_dev*dev0_10m;
dev_0_roth = Ksh*Kth_dev*dev0_10h;

teta_tip = (Dbeta_t12- i_0_rott + dev_0_rott)/(1-m_t+n_LIEt);
teta_mid = (Dbeta_m12- i_0_rotm + dev_0_rotm)/(1-m_m+n_LIEm);
teta_hub = (Dbeta_h12- i_0_roth + dev_0_roth)/(1-m_h+n_LIEh);

i_opt_rott = Ksh*Kth*i0_10tip+n_LIEt*teta_tip;
i_opt_rotm = Ksh*Kth*i0_10mid+n_LIEm*teta_mid;
i_opt_roth = Ksh*Kth*i0_10hub+n_LIEh*teta_hub;

dev_opt_rott = Ksh*Kth_dev*dev0_10t+m_t*teta_tip;
dev_opt_rotm = Ksh*Kth_dev*dev0_10m+m_m*teta_mid;
dev_opt_roth = Ksh*Kth_dev*dev0_10h+m_h*teta_hub;

% CHECK on flow deflection
eps12_t = [Dbeta_t12 teta_tip+i_opt_rott-dev_opt_rott];
eps12_m = [Dbeta_m12 teta_mid+i_opt_rotm-dev_opt_rotm];
eps12_h = [Dbeta_h12 teta_hub+i_opt_roth-dev_opt_roth];

% LIEBLEIN loading criteria
D_t = (w_1t-w_2t)/w_1t+abs((w_1tt-w_2tt)/(2*w_1t*sigma_rt));
D_m = (w_1m-w_2m)/w_1m+abs((w_1tm-w_2tm)/(2*w_1m*sigma_rm));
D_h = (w_1h-w_2h)/w_1h+abs((w_1th-w_2th)/(2*w_1h*sigma_rh));

% The values for D (diffusion factor) proves to be at hub, tip and mid all
% < 0.6 hence there's no need in changing in sigma (?)

% Profile losses according to LIEBLEIN (passing from Johnsen and Bullock) 
YpLJB_t = 0.0035*(1+3.5*D_t+37*(D_t)^4)*2*sigma_rt/cosd(beta_2t);
YpLJB_m = 0.0035*(1+3.5*D_m+37*(D_m)^4)*2*sigma_rm/cosd(beta_2m);
YpLJB_h = 0.0035*(1+3.5*D_h+37*(D_h)^4)*2*sigma_rh/cosd(beta_2h);

% End of the cycle on rotor outlet condition (P_2) and eta_rot
P_tr1t = P_1t*(1+(gamma-1)/2*Mr_1t^2)^(gamma/(gamma-1));
P_tr1m = P_1m*(1+(gamma-1)/2*Mr_1m^2)^(gamma/(gamma-1));
P_tr1h = P_1h*(1+(gamma-1)/2*Mr_1h^2)^(gamma/(gamma-1));

P_tr2t = P_tr1t-YpLJB_t*(P_tr1t-P_1t);
P_tr2m = P_tr1m-YpLJB_m*(P_tr1m-P_1m);
P_tr2h = P_tr1h-YpLJB_h*(P_tr1h-P_1h);

P_2tnew = P_tr2t/((1+(gamma-1)/2*Mr_2t^2)^(gamma/(gamma-1)));
P_2mnew = P_tr2m/((1+(gamma-1)/2*Mr_2m^2)^(gamma/(gamma-1)));
P_2hnew = P_tr2h/((1+(gamma-1)/2*Mr_2h^2)^(gamma/(gamma-1)));

T_2tis = T_1t*(P_2tnew/P_1t)^((gamma-1)/gamma);
T_2mis = T_1m*(P_2mnew/P_1m)^((gamma-1)/gamma);
T_2his = T_1h*(P_2hnew/P_1h)^((gamma-1)/gamma);

eta_rot_tip = [eta_rot_tip (T_2tis-T_1t)/(T_2t-T_1t)];
eta_rot_mid = [eta_rot_mid (T_2mis-T_1m)/(T_2m-T_1m)];
eta_rot_hub = [eta_rot_hub (T_2his-T_1h)/(T_2h-T_1h)];

rho_2t = P_2tnew*10^5/(Rgas/MMa*T_2t);
rho_2m = P_2mnew*10^5/(Rgas/MMa*T_2m);
rho_2h = P_2hnew*10^5/(Rgas/MMa*T_2h);

rho_av_hmt = (rho_2h+rho_2m+rho_2t)/3;
v_2ax = [v_2ax m/(rho_av_hmt*pi*Dmid*b_2)];

w_2ax = v_2ax(end);

o = o+1;
f = f+1;
end


%------------------------------------------------------------------------
%                              STATOR
%------------------------------------------------------------------------

P_t2t = P_2t*(1+(gamma-1)/2*M_2t^2)^(gamma/(gamma-1));
P_t2m = P_2m*(1+(gamma-1)/2*M_2m^2)^(gamma/(gamma-1));
P_t2h = P_2h*(1+(gamma-1)/2*M_2h^2)^(gamma/(gamma-1));

T_t2t = T_2t+v_2t^2/(2*cpa);
T_t2m = T_2m+v_2m^2/(2*cpa);
T_t2h = T_2h+v_2h^2/(2*cpa);

% Notice that the values for T_t2 should be the same along the span since
% according to the free vortex method the work extraction is uniform for
% each radius

T_t3t = T_t2t;
T_t3m = T_t2m;
T_t3h = T_t2h;

%~~~~~ Assumption ~~~~~
% - Repeated stage at midspan v_1m = v_3m; alpha_1m = alpha_3m, etc.
% - Dmid = const.
%~~~~~~~~~~~~~~~~~~~~~~

v_3m = v_1m;
v_3tm = v_1tm;
v_3ax = v_1ax;

alpha_3m = alpha_1m;

P_t3t = [0 P_t3_av]; % first guess
P_t3m = [0 P_t3_av]; % first guess
P_t3h = [0 P_t3_av]; % first guess

tol6 = 1e-4;
e = 1;
while abs(P_t3t(e+1)-P_t3t(e))>tol6 || abs(P_t3m(e+1)-P_t3m(e))>tol6 || abs(P_t3h(e+1)-P_t3h(e))>tol6

T_3m = T_t3m-v_3m^2/(2*cpa);

M_3m = v_3m/sqrt(gamma*Rgas/MMa*T_3m);

P_3m = P_t3m(e+1)/((1+(gamma-1)/2*M_3m^2)^(gamma/(gamma-1)));
rho_3m = P_3m*10^5/(Rgas/MMa*T_3m);

b_3 = m/(rho_3m*pi*Dmid*v_3ax); % different from b_2

Dtip3 = Dmid+b_3/2;
Dhub3 = Dmid-b_3/2;

% Free vortex criterion
% section (3): HUB

v_3th = v_3tm*Dmid/Dhub3;
v_3h = sqrt(v_3th^2+v_3ax^2);

alpha_3h = atand(v_3th/v_3ax);

T_3h = T_t3h-v_3h^2/(2*cpa);

M_3h = v_3h/sqrt(gamma*Rgas/MMa*T_3h);

P_3h = P_t3h(e+1)/((1+(gamma-1)/2*M_3h^2)^(gamma/(gamma-1)));
rho_3h = P_3h*10^5/(Rgas/MMa*T_3h);

% section (3): TIP

v_3tt = v_3tm*Dmid/Dtip3;
v_3t = sqrt(v_3tt^2+v_3ax^2);

alpha_3t = atand(v_3tt/v_3ax);

T_3t = T_t3t-v_3t^2/(2*cpa);

M_3t = v_3t/sqrt(gamma*Rgas/MMa*T_3t);

P_3t = P_t3t(e+1)/((1+(gamma-1)/2*M_3t^2)^(gamma/(gamma-1)));
rho_3t = P_3t*10^5/(Rgas/MMa*T_3t);

% Check on HOWELL correlation: check performed in a graphical way (more
% accurate correlation will be used further on)
Dalpha_t23 = abs(alpha_3t-alpha_2t);
Dalpha_m23 = abs(alpha_3m-alpha_2m);
Dalpha_h23 = abs(alpha_3h-alpha_2h);

% Evaluation of the minimum chord which allows to have a Re = Re_refHOWE =
% 3*10^5

v_av23t = (v_2t+v_3t)/2;
v_av23m = (v_2m+v_3m)/2;
v_av23h = (v_2h+v_3h)/2;

rho_av23t = (rho_2t+rho_3t)/2;
rho_av23m = (rho_2m+rho_3m)/2;
rho_av23h = (rho_2h+rho_3h)/2;

chord_mint_stat = Re_refHOWE*mua/(rho_av23t*v_av23t);
chord_minm_stat = Re_refHOWE*mua/(rho_av23m*v_av23m);
chord_minh_stat = Re_refHOWE*mua/(rho_av23h*v_av23h);

chord_stat = 0.1; % assumed a sthe rotor one

sigma_sm = 1; % (at mid) assumed considering the results previously obtained with HOWELL we have noticed that the blade at midspan is already falling ont the oprimum Howell curve. 
% Moreover this choice for sigma allows to obtain better results even at
% hub and tip since: sigma_sh > 1 reduce the load at the hub (desired since
% overloaded), sigma_st < 1 increase the load at the tip (desired since
% underloaded)

s_statm = chord_stat/sigma_sm;

Nblade_stat = ceil(pi*Dmid/s_statm);

% Since the radius varies along the stator also the pitch, hence the
% solidity varies: we'll consider an average diameter for hub and tip
Dtip_av23 = (Dtip+Dtip3)/2;
Dhub_av23 = (Dhub+Dhub3)/2;

sigma_st = chord_stat/(pi*Dtip_av23/Nblade_stat);
sigma_sh = chord_stat/(pi*Dhub_av23/Nblade_stat);

% LIEBLEIN 

th_c = 0.08; % percentage thickness (assumed equal to rotor)

% Incidence: i_opt_stat = Ksh*Kth*i0_10+n_LIE*teta
Ksh = 1; % for NACA-65

q = 0.28/(0.1+(th_c)^0.3);
Kth = (10*th_c)^q;

p_tip = 0.914+sigma_st^3/160;
p_mid = 0.914+sigma_sm^3/160;
p_hub = 0.914+sigma_sh^3/160;

i0_10tip = abs(alpha_2t)^p_tip/(5+46*exp(-2.3*sigma_st))-0.1*sigma_st^3*exp((abs(alpha_2t)-70)/4);
i0_10mid = abs(alpha_2m)^p_mid/(5+46*exp(-2.3*sigma_sm))-0.1*sigma_sm^3*exp((abs(alpha_2m)-70)/4);
i0_10hub = abs(alpha_2h)^p_hub/(5+46*exp(-2.3*sigma_sh))-0.1*sigma_sh^3*exp((abs(alpha_2h)-70)/4);

n_LIEt = 0.025*sigma_st-0.06-(abs(alpha_2t)/90)^(1+1.2*sigma_st)/(1.5+0.43*sigma_st);
n_LIEm = 0.025*sigma_sm-0.06-(abs(alpha_2m)/90)^(1+1.2*sigma_sm)/(1.5+0.43*sigma_sm);
n_LIEh = 0.025*sigma_sh-0.06-(abs(alpha_2h)/90)^(1+1.2*sigma_sh)/(1.5+0.43*sigma_sh);

i_0_statt = Ksh*Kth*i0_10tip;
i_0_statm = Ksh*Kth*i0_10mid;
i_0_stath = Ksh*Kth*i0_10hub;

% Deviation: dev_opt_stat = Ksh*Kth*dev0_10+m_LIE*teta

x_LIEt = abs(alpha_2t)/100;
x_LIEm = abs(alpha_2m)/100;
x_LIEh = abs(alpha_2h)/100;

m10_t = 0.17-0.0333*x_LIEt+0.333*x_LIEt^2;
m10_m = 0.17-0.0333*x_LIEm+0.333*x_LIEm^2;
m10_h = 0.17-0.0333*x_LIEh+0.333*x_LIEh^2;

b_t = 0.9625-0.17*x_LIEt-0.85*x_LIEt^3;
b_m = 0.9625-0.17*x_LIEm-0.85*x_LIEm^3;
b_h = 0.9625-0.17*x_LIEh-0.85*x_LIEh^3;

Kth_dev = 6.25*(th_c)+37.5*th_c^2;

m_t = m10_t/(sigma_st)^b_t;
m_m = m10_m/(sigma_sm)^b_m;
m_h = m10_h/(sigma_sh)^b_h;

dev0_10t = 0.01*sigma_st*abs(alpha_2t)+(0.74*sigma_st^1.9+3*sigma_st)*(abs(alpha_2t)/90)^(1.67+1.09*sigma_st);
dev0_10m = 0.01*sigma_sm*abs(alpha_2m)+(0.74*sigma_sm^1.9+3*sigma_sm)*(abs(alpha_2m)/90)^(1.67+1.09*sigma_sm);
dev0_10h = 0.01*sigma_sh*abs(alpha_2h)+(0.74*sigma_sh^1.9+3*sigma_sh)*(abs(alpha_2h)/90)^(1.67+1.09*sigma_sh);

dev_0_statt = Ksh*Kth_dev*dev0_10t;
dev_0_statm = Ksh*Kth_dev*dev0_10m;
dev_0_stath = Ksh*Kth_dev*dev0_10h;

teta_tip_stat = (Dalpha_t23- i_0_statt + dev_0_statt)/(1-m_t+n_LIEt);
teta_mid_stat = (Dalpha_m23- i_0_statm + dev_0_statm)/(1-m_m+n_LIEm);
teta_hub_stat = (Dalpha_h23- i_0_stath + dev_0_stath)/(1-m_h+n_LIEh);

i_opt_statt = Ksh*Kth*i0_10tip+n_LIEt*teta_tip_stat;
i_opt_statm = Ksh*Kth*i0_10mid+n_LIEm*teta_mid_stat;
i_opt_stath = Ksh*Kth*i0_10hub+n_LIEh*teta_hub_stat;

dev_opt_statt = Ksh*Kth_dev*dev0_10t+m_t*teta_tip_stat;
dev_opt_statm = Ksh*Kth_dev*dev0_10m+m_m*teta_mid_stat;
dev_opt_stath = Ksh*Kth_dev*dev0_10h+m_h*teta_hub_stat;

% CHECK on flow deflection
eps23_t = [Dalpha_t23 teta_tip_stat+i_opt_statt-dev_opt_statt];
eps23_m = [Dalpha_m23 teta_mid_stat+i_opt_statm-dev_opt_statm];
eps23_h = [Dalpha_h23 teta_hub_stat+i_opt_stath-dev_opt_stath];

% LIEBLEIN loading criteria
D_t_stat = (v_2t-v_3t)/v_2t+abs((v_2tt-v_3tt)/(2*v_2t*sigma_st));
D_m_stat = (v_2m-v_3m)/v_2m+abs((v_2tm-v_3tm)/(2*v_2m*sigma_sm));
D_h_stat = (v_2h-v_3h)/v_2h+abs((v_2th-v_3th)/(2*v_2h*sigma_sh));


% Profile losses according to LIEBLEIN (passing from Johnsen and Bullock) 
YpLJB_t_stat = 0.0035*(1+3.5*D_t_stat+37*(D_t_stat)^4)*2*sigma_st/cosd(alpha_3t);
YpLJB_m_stat = 0.0035*(1+3.5*D_m_stat+37*(D_m_stat)^4)*2*sigma_sm/cosd(alpha_3m);
YpLJB_h_stat = 0.0035*(1+3.5*D_h_stat+37*(D_h_stat)^4)*2*sigma_sh/cosd(alpha_3h);

P_t3t = [P_t3t P_t2t-YpLJB_t_stat*(P_t2t-P_2t)];
P_t3m = [P_t3m P_t2t-YpLJB_m_stat*(P_t2m-P_2m)];
P_t3h = [P_t3h P_t2t-YpLJB_h_stat*(P_t2h-P_2h)];

e = e+1;
end


T_3tis = T_2t*(P_3t/P_2t)^((gamma-1)/gamma);
T_3mis = T_2m*(P_3m/P_2m)^((gamma-1)/gamma);
T_3his = T_2h*(P_3h/P_2h)^((gamma-1)/gamma);

% Stator efficiency (eta_stator_SS)
eta_stat_tip = (T_3tis-T_2t)/(T_3t-T_2t);
eta_stat_mid = (T_3mis-T_2m)/(T_3m-T_2m);
eta_stat_hub = (T_3his-T_2h)/(T_3h-T_2h);

%------------------------------------------------------------------------
%------------------------------------------------------------------------

% Closure of the overall iterative process: check on eta_TT (initially
% assumed from Balje) and P_3av, that must be equal to the required by the
% costumer (P_t3av = beta_tt*P_t0).

P_t3_av = (P_t3t(end)*0.3+P_t3m(end)*0.4+P_t3h(end)*0.3); % arithmetical average
P_t1_av = (P_t1t*0.3+P_t1m*0.4+P_t1h*0.3); % arithmetical average

T_t3is_av = T_t1*(P_t3_av/P_t1_av)^((gamma-1)/gamma); % considering the point t3is as the one obtained with P_3t (real)

eta_tt = [eta_tt (T_t3is_av-T_t0)/(T_t3m-T_t0)]; %T_t3m since T_t3 is constant overall the blade span

g = g+1;

 end

%------------------------------------------------------------------------
%                            RESULTS
%------------------------------------------------------------------------

eta_IGV_SS = (eta_IGV_hub+eta_IGV_mid+eta_IGV_tip)/3; % IGV row efficiency evaluated as an arithmetic average among the values obtained for hub/mid/tip

eta_rot_SS = (eta_rot_hub(end)+eta_rot_mid(end)+eta_rot_tip(end))/3; % Rotor row efficiency evaluated as an arithmetic average among the values obtained for hub/mid/tip

eta_stat_SS = (eta_stat_hub+eta_stat_mid+eta_stat_tip)/3; % stator row efficiency evaluated as an arithmetic average among the values obtained for hub/mid/tip


% Clearance losses by Lakshiminarayana has not been considered, but
% actually as can be proven by the following code lines its contribution is
% almost negligible (order of magnitude is 10^-5)
%
% PSI = Dhs_tt/u_tip^2;
% PHI = m/(rho_1m(end)*u_tip*pi*((Dtip^2-Dhub^2)/4));
% eps = 0.0005; [m] % tip clearance
% alpha_AV12 = atand((tand(alpha_1m)+tand(alpha_2m))/2);
% 
% Deta_clear = 0.07*eps/b_1*PSI/cosd(alpha_AV)*(1+10*sqrt((PHI*eps/chord_rot)/(PSI*cosd(alpha_AV))))



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%                BLADE PROFILES IN DESIGN CONDITION
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% ROTOR
mod_gamma_roth = abs(beta_1h)-abs(teta_hub)/2-(i_opt_roth);
gamma_roth =mod_gamma_roth;

mod_gamma_rotm = abs(beta_1m)-abs(teta_mid)/2-(i_opt_rotm);
gamma_rotm =mod_gamma_rotm;

mod_gamma_rott = abs(beta_1t)-abs(teta_tip)/2-(i_opt_rott);
gamma_rott =mod_gamma_rott;

% STATOR
mod_gamma_stath = abs(alpha_2h)-abs(teta_hub_stat)/2+abs(i_opt_stath);
gamma_stath =mod_gamma_stath;

mod_gamma_statm = abs(alpha_2m)-abs(teta_mid_stat)/2+abs(i_opt_statm);
gamma_statm =mod_gamma_statm;

mod_gamma_statt = abs(alpha_2t)-abs(teta_tip_stat)/2+abs(i_opt_statt);
gamma_statt =mod_gamma_statt;

%-------------------------------------------------------------------------
%                           NACA 65

x_perc=[0 1.25 2.5 5 7.5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percenage x/c
y_t_perc=[0 1.124 1.571 2.222 2.709 3.111 3.746 4.218 4.824 5.057 4.87 4.151 3.038 1.847 0.749 0.354 0.15];% half thickness t/c
y_c_perc=[0 0.535 0.93 1.58 2.12 2.585 3.365 3.98 4.86 5.355 5.515 5.355 4.86 3.98 2.585 1.58 0];% camber line y/c
dy_dx_perc=[Inf 0.3477 0.2915 0.2343 0.1999 0.1749 0.1381 0.1103 0.0675 0.0323 0 0.0323 0.0675 0.1103 0.1749 0.2343 -Inf];% derivative

%------------------------------------------------------------------------
% ROTOR

c_r=chord_rot;

x_r = x_perc*c_r;
y_t_r = c_r*y_t_perc*0.8;
x_rAV = 50*c_r;
y_rAV = 0;

% HUB 

cl_rh=teta_hub/25;

% equations

y_rh = c_r*y_c_perc*cl_rh;
dy_dx_rh = dy_dx_perc*cl_rh;

eps_rh = atand(dy_dx_rh);

xss_rh = x_r-y_t_r.*sind(eps_rh);
yss_rh = y_rh+y_t_r.*cosd(eps_rh);

xps_rh = x_r+y_t_r.*sind(eps_rh);
yps_rh = y_rh-y_t_r.*cosd(eps_rh);

% Passage from x-y coordinates to T-ax

T_coord_ss_rh = -(xss_rh*sind(gamma_roth)+yss_rh*cosd(gamma_roth));
AX_coord_ss_rh = xss_rh*cosd(gamma_roth)-yss_rh*sind(gamma_roth);

T_coord_ps_rh = -(xps_rh*sind(gamma_roth)+yps_rh*cosd(gamma_roth));
AX_coord_ps_rh = xps_rh*cosd(gamma_roth)-yps_rh*sind(gamma_roth);

T_coord_AVrh = -(x_rAV*sind(gamma_roth)+y_rAV*cosd(gamma_roth));
AX_coord_AVrh = x_rAV*cosd(gamma_roth)-y_rAV*sind(gamma_roth);


% MID 

cl_rm=teta_mid/25;

% equations

y_rm = c_r*y_c_perc*cl_rm;
dy_dx_rm = dy_dx_perc*cl_rm;

eps_rm = atand(dy_dx_rm);

xss_rm = x_r-y_t_r.*sind(eps_rm);
yss_rm = y_rm+y_t_r.*cosd(eps_rm);

xps_rm = x_r+y_t_r.*sind(eps_rm);
yps_rm = y_rm-y_t_r.*cosd(eps_rm);

% Passage from x-y coordinates to T-ax

T_coord_ss_rm = -(xss_rm*sind(gamma_rotm)+yss_rm*cosd(gamma_rotm));
AX_coord_ss_rm = xss_rm*cosd(gamma_rotm)-yss_rm*sind(gamma_rotm);

T_coord_ps_rm = -(xps_rm*sind(gamma_rotm)+yps_rm*cosd(gamma_rotm));
AX_coord_ps_rm = xps_rm*cosd(gamma_rotm)-yps_rm*sind(gamma_rotm); 

T_coord_AVrm = -(x_rAV*sind(gamma_rotm)+y_rAV*cosd(gamma_rotm));
AX_coord_AVrm = x_rAV*cosd(gamma_rotm)-y_rAV*sind(gamma_rotm);


% TIP 

cl_rt=teta_tip/25;

% equations

y_rt = c_r*y_c_perc*cl_rt;
dy_dx_rt = dy_dx_perc*cl_rt;

eps_rt = atand(dy_dx_rt);

xss_rt = x_r-y_t_r.*sind(eps_rt);
yss_rt = y_rt+y_t_r.*cosd(eps_rt);

xps_rt = x_r+y_t_r.*sind(eps_rt);
yps_rt = y_rt-y_t_r.*cosd(eps_rt);

% Passage from x-y coordinates to T-ax

T_coord_ss_rt = -(xss_rt*sind(gamma_rott)+yss_rt*cosd(gamma_rott));
AX_coord_ss_rt = xss_rt*cosd(gamma_rott)-yss_rt*sind(gamma_rott);

T_coord_ps_rt = -(xps_rt*sind(gamma_rott)+yps_rt*cosd(gamma_rott));
AX_coord_ps_rt = xps_rt*cosd(gamma_rott)-yps_rt*sind(gamma_rott); 

T_coord_AVrt = -(x_rAV*sind(gamma_rott)+y_rAV*cosd(gamma_rott));
AX_coord_AVrt = x_rAV*cosd(gamma_rott)-y_rAV*sind(gamma_rott);

% plot(xss_rh,yss_rh,xps_rh,yps_rh)

% figure;
% plot(T_coord_ss_rh,AX_coord_ss_rh,T_coord_ps_rh,AX_coord_ps_rh)
% hold on
% 
% plot(T_coord_ss_rm,AX_coord_ss_rm,T_coord_ps_rm,AX_coord_ps_rm)
% hold on
% 
% plot(T_coord_ss_rt,AX_coord_ss_rt,T_coord_ps_rt,AX_coord_ps_rt)
% 
% title('Shape of rotor blade NACA 65 series in x-y plane')
% hold off

% axis([-10 12 -2 12])

% Translation of the blade profiles in T-ax plan in order to correctly
% twist the blade without being inconsistent from mechanical point of view

% HUB
DeltaT_AVrh = T_coord_AVrm - T_coord_AVrh;
DeltaAX_AVrh = AX_coord_AVrm - AX_coord_AVrh;

T_coord_ss_rh_trans = -(xss_rh*sind(gamma_roth)+yss_rh*cosd(gamma_roth))+DeltaT_AVrh;
AX_coord_ss_rh_trans = xss_rh*cosd(gamma_roth)-yss_rh*sind(gamma_roth)+DeltaAX_AVrh;

T_coord_ps_rh_trans = -(xps_rh*sind(gamma_roth)+yps_rh*cosd(gamma_roth))+DeltaT_AVrh;
AX_coord_ps_rh_trans = xps_rh*cosd(gamma_roth)-yps_rh*sind(gamma_roth)+DeltaAX_AVrh;

% TIP
DeltaT_AVrt = T_coord_AVrm - T_coord_AVrt;
DeltaAX_AVrt = AX_coord_AVrm - AX_coord_AVrt;

T_coord_ss_rt_trans = -(xss_rt*sind(gamma_rott)+yss_rt*cosd(gamma_rott))+DeltaT_AVrt;
AX_coord_ss_rt_trans = xss_rt*cosd(gamma_rott)-yss_rt*sind(gamma_rott)+DeltaAX_AVrt;

T_coord_ps_rt_trans = -(xps_rt*sind(gamma_rott)+yps_rt*cosd(gamma_rott))+DeltaT_AVrt;
AX_coord_ps_rt_trans = xps_rt*cosd(gamma_rott)-yps_rt*sind(gamma_rott)+DeltaAX_AVrt;

figure;
plot(T_coord_ss_rh_trans,-AX_coord_ss_rh_trans,'green',T_coord_ps_rh_trans,-AX_coord_ps_rh_trans,'green')
hold on

plot(T_coord_ss_rm,-AX_coord_ss_rm,'red',T_coord_ps_rm,-AX_coord_ps_rm,'red')
hold on

plot(T_coord_ss_rt_trans,-AX_coord_ss_rt_trans,'black',T_coord_ps_rt_trans,-AX_coord_ps_rt_trans,'black')

title('Shape of rotor blade NACA 65 series')
hold off

axis([-10 2 -10 2])

pbaspect([1 1 1])

%------------------------------------------------------------------------
% STATOR

c_s=chord_stat;

x_s = x_perc*c_s;
y_t_s = c_s*y_t_perc*0.8;
x_sAV = 50*c_s;
y_sAV = 0;

% HUB 

cl_sh=teta_hub_stat/25;

% equations

y_sh = c_s*y_c_perc*cl_sh;
dy_dx_sh = dy_dx_perc*cl_sh;

eps_sh = atand(dy_dx_sh);

xss_sh = x_s-y_t_s.*sind(eps_sh);
yss_sh = y_sh+y_t_s.*cosd(eps_sh);

xps_sh = x_s+y_t_s.*sind(eps_sh);
yps_sh = y_sh-y_t_s.*cosd(eps_sh);

% Passage from x-y coordinates to T-ax

T_coord_ss_sh = -(xss_sh*sind(gamma_stath)+yss_sh*cosd(gamma_stath));
AX_coord_ss_sh = xss_sh*cosd(gamma_stath)-yss_sh*sind(gamma_stath);

T_coord_ps_sh = -(xps_sh*sind(gamma_stath)+yps_sh*cosd(gamma_stath));
AX_coord_ps_sh = xps_sh*cosd(gamma_stath)-yps_sh*sind(gamma_stath);

T_coord_AVsh = -(x_sAV*sind(gamma_stath)+y_sAV*cosd(gamma_stath));
AX_coord_AVsh = x_sAV*cosd(gamma_stath)-y_sAV*sind(gamma_stath);


% MID 

cl_sm=teta_mid_stat/25;

% equations

y_sm = c_s*y_c_perc*cl_sm;
dy_dx_sm = dy_dx_perc*cl_sm;

eps_sm = atand(dy_dx_sm);

xss_sm = x_s-y_t_s.*sind(eps_sm);
yss_sm = y_sm+y_t_s.*cosd(eps_sm);

xps_sm = x_s+y_t_s.*sind(eps_sm);
yps_sm = y_sm-y_t_s.*cosd(eps_sm);

% Passage from x-y coordinates to T-ax

T_coord_ss_sm = -(xss_sm*sind(gamma_statm)+yss_sm*cosd(gamma_statm));
AX_coord_ss_sm = xss_sm*cosd(gamma_statm)-yss_sm*sind(gamma_statm);

T_coord_ps_sm = -(xps_sm*sind(gamma_statm)+yps_sm*cosd(gamma_statm));
AX_coord_ps_sm = xps_sm*cosd(gamma_statm)-yps_sm*sind(gamma_statm); 

T_coord_AVsm = -(x_sAV*sind(gamma_statm)+y_sAV*cosd(gamma_statm));
AX_coord_AVsm = x_sAV*cosd(gamma_statm)-y_sAV*sind(gamma_statm);


% TIP 

cl_st=teta_tip_stat/25;

% equations

y_st = c_s*y_c_perc*cl_st;
dy_dx_st = dy_dx_perc*cl_st;

eps_st = atand(dy_dx_st);

xss_st = x_s-y_t_s.*sind(eps_st);
yss_st = y_st+y_t_s.*cosd(eps_st);

xps_st = x_s+y_t_s.*sind(eps_st);
yps_st = y_st-y_t_s.*cosd(eps_st);

% Passage from x-y coordinates to T-ax

T_coord_ss_st = -(xss_st*sind(gamma_statt)+yss_st*cosd(gamma_statt));
AX_coord_ss_st = xss_st*cosd(gamma_statt)-yss_st*sind(gamma_statt);

T_coord_ps_st = -(xps_st*sind(gamma_statt)+yps_st*cosd(gamma_statt));
AX_coord_ps_st = xps_st*cosd(gamma_statt)-yps_st*sind(gamma_statt); 

T_coord_AVst = -(x_sAV*sind(gamma_statt)+y_sAV*cosd(gamma_statt));
AX_coord_AVst = x_sAV*cosd(gamma_statt)-y_sAV*sind(gamma_statt);

% plot(xss_rh,yss_rh,xps_rh,yps_rh)

% figure;
% plot(T_coord_ss_sh,AX_coord_ss_sh,T_coord_ps_sh,AX_coord_ps_sh)
% hold on
% 
% plot(T_coord_ss_sm,AX_coord_ss_sm,T_coord_ps_sm,AX_coord_ps_sm)
% hold on
% 
% plot(T_coord_ss_st,AX_coord_ss_st,T_coord_ps_st,AX_coord_ps_st)
% title('Shape of stator blade NACA 65 series in x-y plane')
% hold off
% 
% axis([-10 12 -2 12])

% Translation of the blade profiles in T-ax plan in order to correctly
% twist the blade without being inconsistent from mechanical point of view

% HUB
DeltaT_AVsh = T_coord_AVsm - T_coord_AVsh;
DeltaAX_AVsh = AX_coord_AVsm - AX_coord_AVsh;

T_coord_ss_sh_trans = -(xss_sh*sind(gamma_stath)+yss_sh*cosd(gamma_stath))+DeltaT_AVsh;
AX_coord_ss_sh_trans = xss_sh*cosd(gamma_stath)-yss_sh*sind(gamma_stath)+DeltaAX_AVsh;

T_coord_ps_sh_trans = -(xps_sh*sind(gamma_stath)+yps_sh*cosd(gamma_stath))+DeltaT_AVsh;
AX_coord_ps_sh_trans = xps_sh*cosd(gamma_stath)-yps_sh*sind(gamma_stath)+DeltaAX_AVsh;

% TIP
DeltaT_AVst = T_coord_AVsm - T_coord_AVst;
DeltaAX_AVst = AX_coord_AVsm - AX_coord_AVst;

T_coord_ss_st_trans = -(xss_st*sind(gamma_statt)+yss_st*cosd(gamma_statt))+DeltaT_AVst;
AX_coord_ss_st_trans = xss_st*cosd(gamma_statt)-yss_st*sind(gamma_statt)+DeltaAX_AVst;

T_coord_ps_st_trans = -(xps_st*sind(gamma_statt)+yps_st*cosd(gamma_statt))+DeltaT_AVst;
AX_coord_ps_st_trans = xps_st*cosd(gamma_statt)-yps_st*sind(gamma_statt)+DeltaAX_AVst;

figure;
plot(-T_coord_ss_sh_trans,-AX_coord_ss_sh_trans,'green',-T_coord_ps_sh_trans,-AX_coord_ps_sh_trans,'green')
hold on

plot(-T_coord_ss_sm,-AX_coord_ss_sm,'red',-T_coord_ps_sm,-AX_coord_ps_sm,'red')
hold on

plot(-T_coord_ss_st_trans,-AX_coord_ss_st_trans,'black',-T_coord_ps_st_trans,-AX_coord_ps_st_trans,'black')

title('Shape of stator blade NACA 65 series')
hold off

axis([-2 10 -10 2])

pbaspect([1 1 1])

%-------------------------------------------------------------------------
%                                IGV

teta_IGV_hub = alpha_1h; % considering that deviation has been neglected and axial inlet
teta_IGV_mid = alpha_1m;
teta_IGV_tip = alpha_1t;

mod_gamma_IGVh = abs(alpha_0)-abs(teta_IGV_hub)/2;
gamma_IGVh =mod_gamma_IGVh;

mod_gamma_IGVm = abs(alpha_0)-abs(teta_IGV_mid)/2;
gamma_IGVm =mod_gamma_IGVm;

mod_gamma_IGVt = abs(alpha_0)-abs(teta_IGV_tip)/2;
gamma_IGVt =mod_gamma_IGVt;

% The blade choosen for the IGV section is the typical NACA A4K6, as also
% suggested by AUNGIER
%-------------------------------------------------------------------------
%                              NACA A4K6

x_perc=[0 1.25 2.5 5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percenage x/c
y_t_perc=[0 0.771 1.057 1.462 2.01 2.386 2.656 2.954 2.971 2.723 2.301 1.87 1.438 1.007 0.791 0];% half thickness t/c
y_c_perc=[0 0.792 1.357 2.248 3.531 4.42 5.04 5.71 5.82 5.516 4.891 4.011 2.922 1.642 0.912 0];% camber line y/c
dy_dx_perc=[Inf 0.5034 0.41 0.3131 0.2110 0.1483 0.1023 0.0359 0.0116 0.0478 0.0761 0.099 0.1184 0.1387 0.155 -Inf];% derivative

%------------------------------------------------------------------------
% ROTOR

c_IGV=chord;

x_IGV = x_perc*c_IGV;
y_t_IGV = c_IGV*y_t_perc*0.8;

% HUB 

cl_IGVh=teta_IGV_hub/25;

% equations

y_IGVh = c_IGV*y_c_perc*cl_IGVh;
dy_dx_IGVh = dy_dx_perc*cl_IGVh;

eps_IGVh = atand(dy_dx_IGVh);

xss_IGVh = x_IGV-y_t_IGV.*sind(eps_IGVh);
yss_IGVh = y_IGVh+y_t_IGV.*cosd(eps_IGVh);

xps_IGVh = x_IGV+y_t_IGV.*sind(eps_IGVh);
yps_IGVh = y_IGVh-y_t_IGV.*cosd(eps_IGVh);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVh = -(xss_IGVh*sind(gamma_IGVh)+yss_IGVh*cosd(gamma_IGVh));
AX_coord_ss_IGVh = xss_IGVh*cosd(gamma_IGVh)-yss_IGVh*sind(gamma_IGVh);

T_coord_ps_IGVh = -(xps_IGVh*sind(gamma_IGVh)+yps_IGVh*cosd(gamma_IGVh));
AX_coord_ps_IGVh = xps_IGVh*cosd(gamma_IGVh)-yps_IGVh*sind(gamma_IGVh);

% MID 

cl_IGVm=teta_IGV_mid/25;

% equations

y_IGVm = c_IGV*y_c_perc*cl_IGVm;
dy_dx_IGVm = dy_dx_perc*cl_IGVm;

eps_IGVm = atand(dy_dx_IGVm);

xss_IGVm = x_IGV-y_t_IGV.*sind(eps_IGVm);
yss_IGVm = y_IGVm+y_t_IGV.*cosd(eps_IGVm);

xps_IGVm = x_IGV+y_t_IGV.*sind(eps_IGVm);
yps_IGVm = y_IGVm-y_t_IGV.*cosd(eps_IGVm);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVm = -(xss_IGVm*sind(gamma_IGVm)+yss_IGVm*cosd(gamma_IGVm));
AX_coord_ss_IGVm = xss_IGVm*cosd(gamma_IGVm)-yss_IGVm*sind(gamma_IGVm);

T_coord_ps_IGVm = -(xps_IGVm*sind(gamma_IGVm)+yps_IGVm*cosd(gamma_IGVm));
AX_coord_ps_IGVm = xps_IGVm*cosd(gamma_IGVm)-yps_IGVm*sind(gamma_IGVm);

% TIP 

cl_IGVt=teta_IGV_tip/25;

% equations

y_IGVt = c_IGV*y_c_perc*cl_IGVt;
dy_dx_IGVt = dy_dx_perc*cl_IGVt;

eps_IGVt = atand(dy_dx_IGVt);

xss_IGVt = x_IGV-y_t_IGV.*sind(eps_IGVt);
yss_IGVt = y_IGVt+y_t_IGV.*cosd(eps_IGVt);

xps_IGVt = x_IGV+y_t_IGV.*sind(eps_IGVt);
yps_IGVt = y_IGVt-y_t_IGV.*cosd(eps_IGVt);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVt = -(xss_IGVt*sind(gamma_IGVt)+yss_IGVt*cosd(gamma_IGVt));
AX_coord_ss_IGVt = xss_IGVt*cosd(gamma_IGVt)-yss_IGVt*sind(gamma_IGVt);

T_coord_ps_IGVt = -(xps_IGVt*sind(gamma_IGVt)+yps_IGVt*cosd(gamma_IGVt));
AX_coord_ps_IGVt = xps_IGVt*cosd(gamma_IGVt)-yps_IGVt*sind(gamma_IGVt); 

% plot(xss_IGVh,yss_IGVh,xps_IGVh,yps_IGVh)

figure;
plot(-T_coord_ss_IGVh,-AX_coord_ss_IGVh,'green',-T_coord_ps_IGVh,-AX_coord_ps_IGVh,'green')
hold on

plot(-T_coord_ss_IGVm,-AX_coord_ss_IGVm,'red',-T_coord_ps_IGVm,-AX_coord_ps_IGVm,'red')
hold on

plot(-T_coord_ss_IGVt,-AX_coord_ss_IGVt,'black',-T_coord_ps_IGVt,-AX_coord_ps_IGVt,'black')

title('Shape of rotor blade NACA A4K6 series in x-y plane')

hold off

% axis([-10 12 -2 12])


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%                   AXIAL COMPRESSOR OFF DESIGN
%------------------------------------------------------------------------
%------------------------------------------------------------------------

m_od = 90; % kg/s

%------------------------------------------------------------------------
%                                IGV
%------------------------------------------------------------------------

% Inlet IGV: iterative cycle to evaluate inlet IGV condition

rho_0_od = [0 rho_t0]; % first guess
z = 1;
tol3 = 1e-4;

while abs(rho_0_od(z+1)-rho_0_od(z))>tol3
v_0_od = m_od/(rho_0_od(z+1)*pi*Dmid*b_0); % constant along the span (axial)
T_0_od = T_t0-v_0_od^2/(2*cpa);
M_0_od = v_0_od/sqrt(gamma*Rgas/MMa*T_0_od);
P_0_od = P_t0/((1+(gamma-1)/2*M_0_od^2)^(gamma/(gamma-1)));

rho_0_od = [rho_0_od P_0_od*10^5/(Rgas/MMa*T_0_od)];

z = z+1;
end

% Outlet IGV: the choice here performed consist in keeping at midspan the
% optimum incidence angle on the rotor inlet, this means that the rotation of the IGV during
% off-design regulation must give a constant beta_1 at the outlet (since
% rotor blades don't move, hence beta_1g doesn't change from design to
% off-design condition). 

%------------------------------------------------------------------------
% Rotor inlet (at MIDSPAN)

v_1axm_od = [0 v_0_od]; % first guess
i = 1;

while abs(v_1axm_od(i+1)-v_1axm_od(i))>tol
w_1m_od = v_1axm_od(i+1)/cosd(beta_1m);
w_1tm_od = w_1m_od*sind(beta_1m);
v_1tm_od = w_1tm_od + u_mid;
v_1m_od = sqrt(v_1axm_od(i+1)^2+v_1tm_od^2);

alpha_1m_od = atand(v_1tm_od/v_1axm_od(i+1));
Dalpha_IGV = alpha_1m_od - alpha_1m; % Dalpha_IGV coincides with the new incidence angle at midspan in off-design condition

T_t1 = T_t0; % IGV is a statoric blade (no work exchanged) even in off design
T_1m_od = T_t1-v_1m_od^2/(2*cpa); % K
M_1m_od = v_1m_od/sqrt(gamma*Rgas/MMa*(T_1m_od));

% AINLEY-MATHIESON correlation for profile losses accounting for incidende
% contribution.
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5

i_IGV_od = Dalpha_IGV;
beta_1prime_od = 90-i_IGV_od;
alpha_2prime_od = 90-alpha_1m_od;

XSI = (90-beta_1prime_od)/(90-alpha_2prime_od);
A = 61.8-(1.6-alpha_2prime_od/165)*alpha_2prime_od;
B = 71.9-1.69*alpha_2prime_od;
C = 7.8-(0.28-alpha_2prime_od/320)*alpha_2prime_od;
D = 14.2-(0.16+alpha_2prime_od/160)*alpha_2prime_od;
i_s0 = 20-(XSI+1)/0.11;

if alpha_2prime_od <= 40
    i_sr = i_s0+A-B*XSI^2+C*XSI^3+D*XSI^4;
else
    XSI40 = (90-beta_1prime_od)/(90-40);
    A40 = 61.8-(1.6-40/165)*40;
    B40 = 71.9-1.69*40;
    C40 = 7.8-(0.28-40/320)*40;
    D40 = 14.2-(0.16+40/160)*40;
    i_s040 = 20-(XSI40+1)/0.11;
    i_sr40 = i_s040+A40-B40*XSI40^2+C40*XSI40^3+D40*XSI40^4;
    i_sr = i_s0+abs(i_sr40-i_s0)*abs(55-alpha_2prime_od)/15;
end

s_c_IGVm_od = s_c_min;
Xi_AM = s_c_IGVm_od-0.75; % s/c keeps the same value as in design condition since the bledes don't change in shape and number

if s_c_IGVm_od <= 0.8
    Dis = -38*Xi_AM-53.5*Xi_AM^2-29*Xi_AM^3;
else
    Dis = 2.0374-(s_c_IGVm_od-0.8)*(69.58-(alpha_2prime_od/14.48)^3.1);
end

is_AM = i_sr+Dis;

if i_IGV_od/is_AM < -3
    Kinc = -1.39214-1.90738*(i_IGV_od/is_AM);
elseif i_IGV_od/is_AM >=-3 && i_IGV_od/is_AM <0
    Kinc = 1+0.52*(abs(i_IGV_od/is_AM))^1.7;
elseif i_IGV_od/is_AM >=0 && i_IGV_od/is_AM <1.7
    Kinc = 1+(i_IGV_od/is_AM)^(2.3+0.5*i_IGV_od/is_AM);
else
    Kinc = 6.23-9.8577*(i_IGV_od/is_AM-1.7);
end

% iterative process for v_1axm_od: the following passages are in accord
% with AINLEY-MATHIESON correlation for the losses

C = 0.08*((alpha_2prime_od/30)^2-1);
n_AM = 1+alpha_2prime_od/30;

if alpha_2prime_od <= 27
    A = 0.025+(27-alpha_2prime_od)/530;
else
    A = 0.025+(27-alpha_2prime_od)/3850;
end

if alpha_2prime_od <= 30
    s_c_minm_od = 0.46+alpha_2prime_od/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2prime_od/1640;
    X_AM = s_c_IGVm_od-s_c_minm_od;
    
    Yp1m_od = A+B*X_AM^2+C*X_AM^3; % proffile losses correlation
else 
    s_c_minm_od = 0.614+alpha_2prime_od/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = s_c_IGVm_od-s_c_minm_od;
    
    Yp1m_od = A+B*(abs(X_AM))^n_AM;
end 

Yp1m_od = Yp1m_od*Kinc; 

syms Pt1m_od P1m_od
sol5 = solve(Yp1m_od-(P_t0-Pt1m_od)/(Pt1m_od-P1m_od)==0, Pt1m_od/P1m_od-(1+(gamma-1)/2*M_1m_od^2)^(gamma/(gamma-1))==0,[Pt1m_od,P1m_od]);
Pt1m_od = double(sol5.Pt1m_od); % to convert the solution from symbolic to numeric
P1m_od = double(sol5.P1m_od); % to convert the solution from symbolic to numeric

rho_1m_od = [ 0 P1m_od*10^5/(Rgas/MMa*T_1m_od)];
j = 1;

while abs(rho_1m_od(j+1)-rho_1m_od(j))>tol
Re_m_od = rho_1m_od(j+1)*v_1m_od*chord/mua; % Re at midspan
Yp1m_Re_od = Yp1m_od*(Re_ref/Re_m_od)^0.2; % Re correction for Yp1

alpha_av01_od = atand((tand(alpha_0)+tand(alpha_1m_od))/2); % alpha_0 is a flow angle hence it doesn't change from design conditions (always 0∞)
cL_od = 2*s_c_IGVm_od*(abs(tand(alpha_1m_od)-tand(alpha_0)))*cosd(alpha_av01_od);
Ysecm_od = chord/b_1vavr*(0.0334*cosd(alpha_1m_od)/cosd(alpha_0))*(cL_od/s_c_IGVm_od)^2*(cosd(alpha_1m_od))^2/(cosd(alpha_av01_od))^3;

Ytotm_od = Yp1m_Re_od+Ysecm_od;

syms Pt1m_od P1m_od
sol6 = solve(Ytotm_od-(P_t0-Pt1m_od)/(Pt1m_od-P1m_od)==0, Pt1m_od/P1m_od-(1+(gamma-1)/2*M_1m_od^2)^(gamma/(gamma-1))==0,[Pt1m_od,P1m_od]);
Pt1m_od = double(sol6.Pt1m_od); % to convert the solution from symbolic to numeric
P1m_od = double(sol6.P1m_od); % to convert the solution from symbolic to numeric

rho_1m_od = [rho_1m_od P1m_od*10^5/(Rgas/MMa*T_1m_od)];

j = j+1;

end

v_1axm_od = [v_1axm_od m_od/(rho_1m_od(end)*pi*Dmid*b_1)];

i = i+1;

end

% vthub = v_1tm_od*Dmid/Dhub;
% alphahub = atand(vthub/v_1axm_od(end))
% alpha_1h_od = alpha_1h + Dalpha_IGV
% 
% vttip = v_1tm_od*Dmid/Dtip;
% alphatip = atand(vttip/v_1axm_od(end))
% alpha_1t_od = alpha_1t + Dalpha_IGV

%------------------------------------------------------------------------
% Rotor inlet (at HUB)
% Both at hub and tip alpha_1 is known since the
% rotation performed by the igv blade is known. Morever we can consider
% that deviation angle is 0∞ (Mach number close to 1) and is not influenced
% by the incidence angle variation at IGV inlet.

v_1axh_od = v_1axm_od(end); % first guess
alpha_1h_od = alpha_1h + Dalpha_IGV;

v_1h_od = v_1axh_od/cosd(alpha_1h_od);
v_1th_od = v_1h_od*sind(alpha_1h_od);
w_1th_od = v_1th_od-u_hub ;
w_1h_od = sqrt(w_1th_od^2+v_1axh_od^2);
beta_1h_od = atand(w_1th_od/v_1axh_od);

T_1h_od = T_t1-v_1h_od^2/(2*cpa); % K
M_1h_od = v_1h_od/sqrt(gamma*Rgas/MMa*(T_1h_od));

% AINLEY-MATHIESON correlation for profile losses accounting for incidende
% contribution.
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5

i_IGV_od = Dalpha_IGV;
beta_1prime_od = 90-i_IGV_od;
alpha_2prime_od = 90-alpha_1h_od;

XSI = (90-beta_1prime_od)/(90-alpha_2prime_od);
A = 61.8-(1.6-alpha_2prime_od/165)*alpha_2prime_od;
B = 71.9-1.69*alpha_2prime_od;
C = 7.8-(0.28-alpha_2prime_od/320)*alpha_2prime_od;
D = 14.2-(0.16+alpha_2prime_od/160)*alpha_2prime_od;
i_s0 = 20-(XSI+1)/0.11;

if alpha_2prime_od <= 40
    i_sr = i_s0+A-B*XSI^2+C*XSI^3+D*XSI^4;
else
    XSI40 = (90-beta_1prime_od)/(90-40);
    A40 = 61.8-(1.6-40/165)*40;
    B40 = 71.9-1.69*40;
    C40 = 7.8-(0.28-40/320)*40;
    D40 = 14.2-(0.16+40/160)*40;
    i_s040 = 20-(XSI40+1)/0.11;
    i_sr40 = i_s040+A40-B40*XSI40^2+C40*XSI40^3+D40*XSI40^4;
    i_sr = i_s0+abs(i_sr40-i_s0)*abs(55-alpha_2prime_od)/15;
end

s_c_IGVh_od = s_c_h;
Xi_AM = s_c_IGVh_od-0.75; % s/c keeps the same value as in design condition since the bledes don't change in shape and number

if s_c_IGVh_od <= 0.8
    Dis = -38*Xi_AM-53.5*Xi_AM^2-29*Xi_AM^3;
else
    Dis = 2.0374-(s_c_IGVh_od-0.8)*(69.58-(alpha_2prime_od/14.48)^3.1);
end

is_AM = i_sr+Dis;

if i_IGV_od/is_AM < -3
    Kinc = -1.39214-1.90738*(i_IGV_od/is_AM);
elseif i_IGV_od/is_AM >=-3 && i_IGV_od/is_AM <0
    Kinc = 1+0.52*(abs(i_IGV_od/is_AM))^1.7;
elseif i_IGV_od/is_AM >=0 && i_IGV_od/is_AM <1.7
    Kinc = 1+(i_IGV_od/is_AM)^(2.3+0.5*i_IGV_od/is_AM);
else
    Kinc = 6.23-9.8577*(i_IGV_od/is_AM-1.7);
end

% iterative process for v_1axh_od: the following passages are in accord
% with AINLEY-MATHIESON correlation for the losses

C = 0.08*((alpha_2prime_od/30)^2-1);
n_AM = 1+alpha_2prime_od/30;

if alpha_2prime_od <= 27
    A = 0.025+(27-alpha_2prime_od)/530;
else
    A = 0.025+(27-alpha_2prime_od)/3850;
end

if alpha_2prime_od <= 30
    s_c_minh_od = 0.46+alpha_2prime_od/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2prime_od/1640;
    X_AM = s_c_IGVh_od-s_c_minh_od;
    
    Yp1h_od = A+B*X_AM^2+C*X_AM^3; % proffile losses correlation
else 
    s_c_minh_od = 0.614+alpha_2prime_od/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = s_c_IGVh_od-s_c_minh_od;
    
    Yp1h_od = A+B*(abs(X_AM))^n_AM;
end 

Yp1h_od = Yp1h_od*Kinc; 

syms Pt1h_od P1h_od
sol7 = solve(Yp1h_od-(P_t0-Pt1h_od)/(Pt1h_od-P1h_od)==0, Pt1h_od/P1h_od-(1+(gamma-1)/2*M_1h_od^2)^(gamma/(gamma-1))==0,[Pt1h_od,P1h_od]);
Pt1h_od = double(sol7.Pt1h_od); % to convert the solution from symbolic to numeric
P1h_od = double(sol7.P1h_od); % to convert the solution from symbolic to numeric

rho_1h_od = [0 P1h_od*10^5/(Rgas/MMa*T_1h_od)];
i = 1;

while abs(rho_1h_od(i+1)-rho_1h_od(i))>tol
Re_h_od = rho_1h_od(i+1)*v_1h_od*chord/mua; % Re at midspan
Yp1h_Re_od = Yp1h_od*(Re_ref/Re_h_od)^0.2; % Re correction for Yp1

alpha_av01_od = atand((tand(alpha_0)+tand(alpha_1h_od))/2); % alpha_0 is a flow angle hence it doesn't change from design conditions (always 0∞)
cL_od = 2*s_c_IGVh_od*(abs(tand(alpha_1h_od)-tand(alpha_0)))*cosd(alpha_av01_od);
Ysech_od = chord/b_1vavr*(0.0334*cosd(alpha_1h_od)/cosd(alpha_0))*(cL_od/s_c_IGVh_od)^2*(cosd(alpha_1h_od))^2/(cosd(alpha_av01_od))^3;

Ytoth_od = Yp1h_Re_od+Ysech_od;

syms Pt1h_od P1h_od
sol8 = solve(Ytoth_od-(P_t0-Pt1h_od)/(Pt1h_od-P1h_od)==0, Pt1h_od/P1h_od-(1+(gamma-1)/2*M_1h_od^2)^(gamma/(gamma-1))==0,[Pt1h_od,P1h_od]);
Pt1h_od = double(sol8.Pt1h_od); % to convert the solution from symbolic to numeric
P1h_od = double(sol8.P1h_od); % to convert the solution from symbolic to numeric

rho_1h_od = [rho_1h_od P1h_od*10^5/(Rgas/MMa*T_1h_od)];

i = i+1;

end


%------------------------------------------------------------------------
% Rotor inlet (at TIP)
% Both at hub and tip alpha_1 is known since the
% rotation performed by the igv blade is known. Morever we can consider
% that deviation angle is 0∞ (Mach number close to 1) and is not influenced
% by the incidence angle variation at IGV inlet.

v_1axt_od = v_1axm_od(end); % first guess
alpha_1t_od = alpha_1t + Dalpha_IGV;

v_1t_od = v_1axt_od/cosd(alpha_1t_od);
v_1tt_od = v_1t_od*sind(alpha_1t_od);
w_1tt_od = v_1tt_od-u_tip ;
w_1t_od = sqrt(w_1tt_od^2+v_1axt_od^2);
beta_1t_od = atand(w_1tt_od/v_1axt_od);

T_1t_od = T_t1-v_1t_od^2/(2*cpa); % K
M_1t_od = v_1t_od/sqrt(gamma*Rgas/MMa*(T_1t_od));

% AINLEY-MATHIESON correlation for profile losses accounting for incidende
% contribution.
% alpha_2prime = 90∞-alpha_1
% beta_1prime = 90∞-alpha_0
% Re_ref = 2*10^5

i_IGV_od = Dalpha_IGV;
beta_1prime_od = 90-i_IGV_od;
alpha_2prime_od = 90-alpha_1t_od;

XSI = (90-beta_1prime_od)/(90-alpha_2prime_od);
A = 61.8-(1.6-alpha_2prime_od/165)*alpha_2prime_od;
B = 71.9-1.69*alpha_2prime_od;
C = 7.8-(0.28-alpha_2prime_od/320)*alpha_2prime_od;
D = 14.2-(0.16+alpha_2prime_od/160)*alpha_2prime_od;
i_s0 = 20-(XSI+1)/0.11;

if alpha_2prime_od <= 40
    i_sr = i_s0+A-B*XSI^2+C*XSI^3+D*XSI^4;
else
    XSI40 = (90-beta_1prime_od)/(90-40);
    A40 = 61.8-(1.6-40/165)*40;
    B40 = 71.9-1.69*40;
    C40 = 7.8-(0.28-40/320)*40;
    D40 = 14.2-(0.16+40/160)*40;
    i_s040 = 20-(XSI40+1)/0.11;
    i_sr40 = i_s040+A40-B40*XSI40^2+C40*XSI40^3+D40*XSI40^4;
    i_sr = i_s0+abs(i_sr40-i_s0)*abs(55-alpha_2prime_od)/15;
end

s_c_IGVt_od = s_c_t;
Xi_AM = s_c_IGVt_od-0.75; % s/c keeps the same value as in design condition since the bledes don't change in shape and number

if s_c_IGVt_od <= 0.8
    Dis = -38*Xi_AM-53.5*Xi_AM^2-29*Xi_AM^3;
else
    Dis = 2.0374-(s_c_IGVt_od-0.8)*(69.58-(alpha_2prime_od/14.48)^3.1);
end

is_AM = i_sr+Dis;

if i_IGV_od/is_AM < -3
    Kinc = -1.39214-1.90738*(i_IGV_od/is_AM);
elseif i_IGV_od/is_AM >=-3 && i_IGV_od/is_AM <0
    Kinc = 1+0.52*(abs(i_IGV_od/is_AM))^1.7;
elseif i_IGV_od/is_AM >=0 && i_IGV_od/is_AM <1.7
    Kinc = 1+(i_IGV_od/is_AM)^(2.3+0.5*i_IGV_od/is_AM);
else
    Kinc = 6.23-9.8577*(i_IGV_od/is_AM-1.7);
end

% iterative process for v_1axh_od: the following passages are in accord
% with AINLEY-MATHIESON correlation for the losses

C = 0.08*((alpha_2prime_od/30)^2-1);
n_AM = 1+alpha_2prime_od/30;

if alpha_2prime_od <= 27
    A = 0.025+(27-alpha_2prime_od)/530;
else
    A = 0.025+(27-alpha_2prime_od)/3850;
end

if alpha_2prime_od <= 30
    s_c_mint_od = 0.46+alpha_2prime_od/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2prime_od/1640;
    X_AM = s_c_IGVt_od-s_c_mint_od;
    
    Yp1t_od = A+B*X_AM^2+C*X_AM^3; % proffile losses correlation
else 
    s_c_mint_od = 0.614+alpha_2prime_od/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = s_c_IGVt_od-s_c_mint_od;
    
    Yp1t_od = A+B*(abs(X_AM))^n_AM;
end 

Yp1t_od = Yp1t_od*Kinc; 

syms Pt1t_od P1t_od
sol9 = solve(Yp1t_od-(P_t0-Pt1t_od)/(Pt1t_od-P1t_od)==0, Pt1t_od/P1t_od-(1+(gamma-1)/2*M_1t_od^2)^(gamma/(gamma-1))==0,[Pt1t_od,P1t_od]);
Pt1t_od = double(sol9.Pt1t_od); % to convert the solution from symbolic to numeric
P1t_od = double(sol9.P1t_od); % to convert the solution from symbolic to numeric

rho_1t_od = [0 P1t_od*10^5/(Rgas/MMa*T_1t_od)];
i = 1;

while abs(rho_1t_od(i+1)-rho_1t_od(i))>tol 
Re_t_od = rho_1t_od(i+1)*v_1t_od*chord/mua; % Re at midspan
Yp1t_Re_od = Yp1t_od*(Re_ref/Re_t_od)^0.2; % Re correction for Yp1

alpha_av01_od = atand((tand(alpha_0)+tand(alpha_1t_od))/2); % alpha_0 is a flow angle hence it doesn't change from design conditions (always 0∞)
cL_od = 2*s_c_IGVt_od*(abs(tand(alpha_1t_od)-tand(alpha_0)))*cosd(alpha_av01_od);
Ysect_od = chord/b_1vavr*(0.0334*cosd(alpha_1t_od)/cosd(alpha_0))*(cL_od/s_c_IGVt_od)^2*(cosd(alpha_1t_od))^2/(cosd(alpha_av01_od))^3;

Ytott_od = Yp1t_Re_od+Ysect_od;

syms Pt1t_od P1t_od
sol10 = solve(Ytott_od-(P_t0-Pt1t_od)/(Pt1t_od-P1t_od)==0, Pt1t_od/P1t_od-(1+(gamma-1)/2*M_1t_od^2)^(gamma/(gamma-1))==0,[Pt1t_od,P1t_od]);
Pt1t_od = double(sol10.Pt1t_od); % to convert the solution from symbolic to numeric
P1t_od = double(sol10.P1t_od); % to convert the solution from symbolic to numeric

rho_1t_od = [rho_1t_od P1t_od*10^5/(Rgas/MMa*T_1t_od)];

i = i+1;

end

Mr_1m_od = w_1m_od/sqrt(Rgas/MMa*gamma*T_1m_od);

%------------------------------------------------------------------------
%                          ROTOR OUTLET
%------------------------------------------------------------------------

v_2axm_od = [0 v_1axm_od(end)]; % first guess
j = 1;

while abs(v_2axm_od(j+1)-v_2axm_od(j))>tol
w_2m_od = v_2axm_od(j+1)/cosd(beta_2m); % considering JOHNSEN AND BULLOCK correlation for 
% the variation of deviation angle in relation to the change of
% incidence angle, we can consider that if i = i_opt (as in the case of
% midspan), dev = dev_opt, hence, since beta_2geom doesn't change also
% beta_2m keeps constant.

w_2tm_od = w_2m_od*sind(beta_2m);
v_2tm_od = w_2tm_od+u_mid;
alpha_2m_od = atand(v_2tm_od/v_2axm_od(j+1));
v_2m_od = sqrt(v_2tm_od^2+v_2axm_od(j+1)^2);

Leul_od = u_mid*(v_2tm_od-v_1tm_od);
T_t2_od = T_t1+Leul_od/cpa;
T_2m_od = T_t2_od-v_2m_od^2/(2*cpa);

% LIEBLEIN loading criteria
D_m_od = (w_1m_od-w_2m_od)/w_1m_od+abs((w_1tm_od-w_2tm_od)/(2*w_1m_od*sigma_rm));
YpLJB_m_od = 0.0035*(1+3.5*D_m_od+37*(D_m_od)^4)*2*sigma_rm/cosd(beta_2m);

P_tr1m_od = P1m_od*(1+(gamma-1)/2*Mr_1m_od^2)^(gamma/(gamma-1));
P_tr2m_od = P_tr1m_od-YpLJB_m_od*(P_tr1m_od-P1m_od);

Mr_2m_od = w_2m_od/sqrt(Rgas/MMa*gamma*T_2m_od);
P2m_od= P_tr2m_od/(1+(gamma-1)/2*Mr_2m_od^2)^(gamma/(gamma-1));
rho_2m_od = P2m_od*10^5/(Rgas/MMa*T_2m_od);

v_2axm_od = [v_2axm_od m_od/(rho_2m_od*pi*Dmid*b_2)];

j=j+1;
end

M_2m_od = v_2m_od/sqrt(gamma*Rgas/MMa*T_2m_od);

%------------------------------------------------------------------------
% Rotor outlet (at tip)

% Losses correlation proposed by AUNGIER
e = 0.65-0.002*teta_tip;
alpha_star =(3.6*Ksh*Kth+0.3532*teta_tip*(0.5)^0.25)*sigma_rt^e; % Correlation proposed by AUNGIER (a/c=0.5 for NACA-65)

syms alpha_c
sol10 = solve(alpha_c-alpha_star-(-9+(1-(30/(alpha_c+gamma_rott)^0.48)*teta_tip/4.176))==0,[alpha_c]);
alpha_c = max(real(double(sol10))); % to convert the solution from symbolic to numeric
% Please notice that the right solution is the maximum one since AUNGIER
% suggest to take beta_1c = alpha_c+gamma_rott hiher than 20∞

syms alpha_s
sol10 = solve(alpha_s-alpha_star-(10.3+(2.92-(alpha_s+gamma_rott)/15.6)*teta_tip/8.2)==0,[alpha_s]);
alpha_s = double(sol10);

Rc = alpha_star-alpha_c;
Rs = alpha_s-alpha_star;

Mr_1t_od = w_1t_od/sqrt(Rgas/MMa*gamma*T_1t_od);

i_c = i_opt_rott-Rc/(1+0.5*Mr_1t_od^3);
i_s = i_opt_rott+Rs/(1+0.5*(Ksh*Mr_1t_od)^3);
i_m = i_c+(i_s-i_c)*Rc/(Rc+Rs);

omega_m = YpLJB_t*(1+(i_m-i_opt_rott)^2/Rs^2);

i_tipr_od = beta_1t_od-(beta_1t-i_opt_rott);

if i_tipr_od >= i_m
    XSI_i = (i_tipr_od-i_m)/(i_s-i_m);
else 
    XSI_i = (i_tipr_od-i_m)/(i_m-i_c);
end
    
if XSI_i >= -2 && XSI_i <= 1
    YpLJB_t_od = omega_m*(1+XSI_i^2);
elseif XSI <-2
    YpLJB_t_od = omega_m*(5-4*(XSI_i+2));
else
    YpLJB_t_od = omega_m*(2+2*(XSI_i-1));
end

%free vortex
v_2tt_od = v_2tm_od*Dmid/Dtip;
v_2t_od = sqrt(v_2tt_od^2+v_2axm_od(end)^2);
w_2tt_od = v_2tt_od-u_tip;
w_2t_od = sqrt(w_2tt_od^2+v_2axm_od(end)^2);

alpha_2t_od = atand(v_2tt_od/v_2axm_od(end));
beta_2t_od = atand(w_2tt_od/v_2axm_od(end));    
    
T_2t_od = T_t2_od-v_2t_od^2/(2*cpa);

P_tr1t_od = P1t_od*(1+(gamma-1)/2*Mr_1t_od^2)^(gamma/(gamma-1));
P_tr2t_od = P_tr1t_od-YpLJB_t_od*(P_tr1t_od-P1t_od);

Mr_2t_od = w_2t_od/sqrt(Rgas/MMa*gamma*T_2t_od);
P2t_od= P_tr2t_od/(1+(gamma-1)/2*Mr_2t_od^2)^(gamma/(gamma-1));
rho_2t_od = P2t_od*10^5/(Rgas/MMa*T_2t_od);

%------------------------------------------------------------------------
% Rotor outlet (at hub)

% Losses correlation proposed by AUNGIER
e = 0.65-0.002*teta_hub;
alpha_starh =(3.6*Ksh*Kth+0.3532*teta_hub*(0.5)^0.25)*sigma_rh^e; % Correlation proposed by AUNGIER (a/c=0.5 for NACA-65)

syms alpha_c
sol11 = solve(alpha_c-alpha_star-(-9+(1-(30/(alpha_c+gamma_roth)^0.48)*teta_hub/4.176))==0,[alpha_c]);
alpha_c = max(real(double(sol11))); % to convert the solution from symbolic to numeric
% Please notice that the right solution is the maximum one since AUNGIER
% suggest to take beta_1c = alpha_c+gamma_rott hiher than 20∞

syms alpha_s
sol12 = solve(alpha_s-alpha_star-(10.3+(2.92-(alpha_s+gamma_roth)/15.6)*teta_hub/8.2)==0,[alpha_s]);
alpha_s = double(sol12);

Rc = alpha_star-alpha_c;
Rs = alpha_s-alpha_star;

Mr_1h_od = w_1h_od/sqrt(Rgas/MMa*gamma*T_1h_od);

i_c = i_opt_roth-Rc/(1+0.5*Mr_1h_od^3);
i_s = i_opt_roth+Rs/(1+0.5*(Ksh*Mr_1h_od)^3);
i_m = i_c+(i_s-i_c)*Rc/(Rc+Rs);

omega_m = YpLJB_h*(1+(i_m-i_opt_roth)^2/Rs^2);

i_hubr_od = beta_1h_od-(beta_1h-i_opt_roth);

if i_hubr_od >= i_m
    XSI_i = (i_hubr_od-i_m)/(i_s-i_m);
else 
    XSI_i = (i_hubr_od-i_m)/(i_m-i_c);
end
    
if XSI_i >= -2 && XSI_i <= 1
    YpLJB_h_od = omega_m*(1+XSI_i^2);
elseif XSI <-2
    YpLJB_h_od = omega_m*(5-4*(XSI_i+2));
else
    YpLJB_h_od = omega_m*(2+2*(XSI_i-1));
end

%free vortex
v_2th_od = v_2tm_od*Dmid/Dhub;
v_2h_od = sqrt(v_2th_od^2+v_2axm_od(end)^2);
w_2th_od = v_2th_od-u_hub;
w_2h_od = sqrt(w_2th_od^2+v_2axm_od(end)^2);

alpha_2h_od = atand(v_2th_od/v_2axm_od(end));
beta_2h_od = atand(w_2th_od/v_2axm_od(end));    
    
T_2h_od = T_t2_od-v_2h_od^2/(2*cpa);

P_tr1h_od = P1h_od*(1+(gamma-1)/2*Mr_1h_od^2)^(gamma/(gamma-1));
P_tr2h_od = P_tr1h_od-YpLJB_h_od*(P_tr1h_od-P1h_od);

Mr_2h_od = w_2h_od/sqrt(Rgas/MMa*gamma*T_2h_od);
P2h_od= P_tr2h_od/(1+(gamma-1)/2*Mr_2h_od^2)^(gamma/(gamma-1));
rho_2h_od = P2h_od*10^5/(Rgas/MMa*T_2h_od);
    

%------------------------------------------------------------------------
%                          STATOR OUTLET
%------------------------------------------------------------------------   

% Stator outlet (at MID)

T_t3_od = T_t2_od; % Since statoric component
v_3axm_od = [0 v_2axm_od(end)]; % first guess
i = 1;

% Considering now the correlation providing the change of deviation angle
% in accordance to the variation of incidence angle provided by JOHNSEN AND
% BULLOK

while abs(v_3axm_od(i+1)-v_3axm_od(i))>tol
diff_dev_i = (1+(sigma_sm+0.25*sigma_sm)^4*(alpha_2m_od/53)^2.5)/exp(3.1*sigma_sm);
i_mids_od = alpha_2m_od-(alpha_2m-i_opt_statm);
dev_mids_od = dev_opt_statm+diff_dev_i*(i_mids_od-i_opt_statm)+10*(1-v_3axm_od(i+1)/v_2axm_od(end));

alpha_3m_od = alpha_3m-dev_opt_statm+dev_mids_od;

v_3m_od = v_3axm_od(i+1)/cosd(alpha_3m_od);
v_3tm_od = v_3m_od*sind(alpha_3m_od);

T_3m_od = T_t3_od-v_3m_od^2/(2*cpa);
M_3m_od = v_3m_od/sqrt(gamma*Rgas/MMa*T_3m_od);

% Losses correlation proposed by AUNGIER
e = 0.65-0.002*teta_mid_stat;
alpha_star =(3.6*Ksh*Kth+0.3532*teta_mid_stat*(0.5)^0.25)*sigma_sm^e; % Correlation proposed by AUNGIER (a/c=0.5 for NACA-65)

syms alpha_c
sol13 = solve(alpha_c-alpha_star-(-9+(1-(30/(alpha_c+gamma_statm)^0.48)*teta_mid_stat/4.176))==0,[alpha_c]);
alpha_c = max(real(double(sol13))); % to convert the solution from symbolic to numeric
% Please notice that the right solution is the maximum one since AUNGIER
% suggest to take beta_1c = alpha_c+gamma_rott hiher than 20∞

syms alpha_s
sol14 = solve(alpha_s-alpha_star-(10.3+(2.92-(alpha_s+gamma_statm)/15.6)*teta_mid_stat/8.2)==0,[alpha_s]);
alpha_s = double(sol14);

Rc = alpha_star-alpha_c;
Rs = alpha_s-alpha_star;

i_c = i_opt_statm-Rc/(1+0.5*M_2m_od^3);
i_s = i_opt_statm+Rs/(1+0.5*(Ksh*M_2m_od)^3);
i_m = i_c+(i_s-i_c)*Rc/(Rc+Rs);

omega_m = YpLJB_m_stat*(1+(i_m-i_opt_statm)^2/Rs^2);

if i_mids_od >= i_m
    XSI_i = (i_mids_od-i_m)/(i_s-i_m);
else 
    XSI_i = (i_mids_od-i_m)/(i_m-i_c);
end
    
if XSI_i >= -2 && XSI_i <= 1
    YpLJB_m_statod = omega_m*(1+XSI_i^2);
elseif XSI <-2
    YpLJB_m_statod = omega_m*(5-4*(XSI_i+2));
else
    YpLJB_m_statod = omega_m*(2+2*(XSI_i-1));
end  

P_t2m_od = P2m_od*(1+(gamma-1)/2*M_2m_od^2)^(gamma/(gamma-1));
P_t3m_od = P_t2m_od-YpLJB_m_statod*(P_t2m_od-P2m_od);

P3m_od= P_t3m_od/(1+(gamma-1)/2*M_3m_od^2)^(gamma/(gamma-1));
rho_3m_od = P3m_od*10^5/(Rgas/MMa*T_3m_od);

v_3axm_od = [v_3axm_od m_od/(rho_3m_od*pi*Dmid*b_3)];

i = i+1;

end

%------------------------------------------------------------------------
% Stator outlet (at TIP)

% Losses correlation proposed by AUNGIER
e = 0.65-0.002*teta_tip_stat;
alpha_star =(3.6*Ksh*Kth+0.3532*teta_tip_stat*(0.5)^0.25)*sigma_st^e; % Correlation proposed by AUNGIER (a/c=0.5 for NACA-65)

syms alpha_c
sol15 = solve(alpha_c-alpha_star-(-9+(1-(30/(alpha_c+gamma_statt)^0.48)*teta_tip_stat/4.176))==0,[alpha_c]);
alpha_c = max(real((double(sol15)))); % to convert the solution from symbolic to numeric
% Please notice that the right solution is the maximum one since AUNGIER
% suggest to take beta_1c = alpha_c+gamma_rott hiher than 20∞

syms alpha_s
sol16 = solve(alpha_s-alpha_star-(10.3+(2.92-(alpha_s+gamma_statt)/15.6)*teta_tip_stat/8.2)==0,[alpha_s]);
alpha_s = double(sol16);

Rc = alpha_star-alpha_c;
Rs = alpha_s-alpha_star;

M_2t_od = v_2t_od/sqrt(Rgas/MMa*gamma*T_2t_od);

i_c = i_opt_statt-Rc/(1+0.5*M_2t_od^3);
i_s = i_opt_statt+Rs/(1+0.5*(Ksh*M_2t_od)^3);
i_m = i_c+(i_s-i_c)*Rc/(Rc+Rs);

omega_m = YpLJB_t_stat*(1+(i_m-i_opt_statt)^2/Rs^2);

i_tips_od = alpha_2t_od-(alpha_2t-i_opt_statt);

if i_tips_od >= i_m
    XSI_i = (i_tips_od-i_m)/(i_s-i_m);
else 
    XSI_i = (i_tips_od-i_m)/(i_m-i_c);
end
    
if XSI_i >= -2 && XSI_i <= 1
    YpLJB_t_statod = omega_m*(1+XSI_i^2);
elseif XSI <-2
    YpLJB_t_statod = omega_m*(5-4*(XSI_i+2));
else
    YpLJB_t_statod = omega_m*(2+2*(XSI_i-1));
end

%free vortex
v_3tt_od = v_3tm_od*Dmid/Dtip3;
v_3t_od = sqrt(v_3tt_od^2+v_3axm_od(end)^2);

alpha_3t_od = atand(v_3tt_od/v_3axm_od(end));

    
T_3t_od = T_t3_od-v_3t_od^2/(2*cpa);
M_3t_od = v_3t_od/sqrt(gamma*Rgas/MMa*T_3t_od);

P_t2t_od = P2t_od*(1+(gamma-1)/2*M_2t_od^2)^(gamma/(gamma-1));
P_t3t_od = P_t2t_od-YpLJB_t_statod*(P_t2t_od-P2t_od);

P3t_od= P_t3t_od/(1+(gamma-1)/2*M_3t_od^2)^(gamma/(gamma-1));
rho_3t_od = P3t_od*10^5/(Rgas/MMa*T_3t_od);

%------------------------------------------------------------------------
% Stator outlet (at HUB)

% Losses correlation proposed by AUNGIER
e = 0.65-0.002*teta_hub_stat;
alpha_star =(3.6*Ksh*Kth+0.3532*teta_hub_stat*(0.5)^0.25)*sigma_sh^e; % Correlation proposed by AUNGIER (a/c=0.5 for NACA-65)

syms alpha_c
sol16 = solve(alpha_c-alpha_star-(-9+(1-(30/(alpha_c+gamma_stath)^0.48)*teta_hub_stat/4.176))==0,[alpha_c]);
alpha_c = max(real((double(sol16)))); % to convert the solution from symbolic to numeric
% Please notice that the right solution is the maximum one since AUNGIER
% suggest to take beta_1c = alpha_c+gamma_rott hiher than 20∞

syms alpha_s
sol17 = solve(alpha_s-alpha_star-(10.3+(2.92-(alpha_s+gamma_stath)/15.6)*teta_hub_stat/8.2)==0,[alpha_s]);
alpha_s = double(sol17);

Rc = alpha_star-alpha_c;
Rs = alpha_s-alpha_star;

M_2h_od = v_2h_od/sqrt(Rgas/MMa*gamma*T_2h_od);

i_c = i_opt_stath-Rc/(1+0.5*M_2h_od^3);
i_s = i_opt_stath+Rs/(1+0.5*(Ksh*M_2h_od)^3);
i_m = i_c+(i_s-i_c)*Rc/(Rc+Rs);

omega_m = YpLJB_h_stat*(1+(i_m-i_opt_stath)^2/Rs^2);

i_hubs_od = alpha_2h_od-(alpha_2h-i_opt_stath);

if i_hubs_od >= i_m
    XSI_i = (i_hubs_od-i_m)/(i_s-i_m);
else 
    XSI_i = (i_hubs_od-i_m)/(i_m-i_c);
end
    
if XSI_i >= -2 && XSI_i <= 1
    YpLJB_h_statod = omega_m*(1+XSI_i^2);
elseif XSI <-2
    YpLJB_h_statod = omega_m*(5-4*(XSI_i+2));
else
    YpLJB_h_statod = omega_m*(2+2*(XSI_i-1));
end

%free vortex
v_3th_od = v_3tm_od*Dmid/Dhub3;
v_3h_od = sqrt(v_3th_od^2+v_3axm_od(end)^2);

alpha_3h_od = atand(v_3th_od/v_3axm_od(end));

    
T_3h_od = T_t3_od-v_3h_od^2/(2*cpa);
M_3h_od = v_3h_od/sqrt(gamma*Rgas/MMa*T_3h_od);

P_t2h_od = P2h_od*(1+(gamma-1)/2*M_2h_od^2)^(gamma/(gamma-1));
P_t3h_od = P_t2h_od-YpLJB_h_statod*(P_t2h_od-P2h_od);

P3h_od= P_t3h_od/(1+(gamma-1)/2*M_3h_od^2)^(gamma/(gamma-1));
rho_3h_od = P3h_od*10^5/(Rgas/MMa*T_3h_od);

%------------------------------------------------------------------------
%                            RESULTS
%------------------------------------------------------------------------

beta_tt_od = P_t3m_od/P_t0;
Dhs_tt_od = cpa*T_t0*(beta_tt_od^((gamma-1)/gamma)-1);
eta_tt_od = Dhs_tt_od/Leul_od;


%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%                                IGV OFF-DESIGN
%------------------------------------------------------------------------
figure;
plot(-T_coord_ss_IGVm,-AX_coord_ss_IGVm,'black',-T_coord_ps_IGVm,-AX_coord_ps_IGVm,'black')
hold on
% This figure will be completed at the end of this section and will allow
% to study the position of blade profile at midspan in off-design condition
% with respect to the one in design.
%-------------------------------------------------------------------------

teta_IGV_hub = alpha_1h; % considering that deviation has been neglected and axial inlet
teta_IGV_mid = alpha_1m;
teta_IGV_tip = alpha_1t;

mod_gamma_IGVh = abs(alpha_0)-abs(teta_IGV_hub)/2-i_IGV_od;
gamma_IGVh =mod_gamma_IGVh;

mod_gamma_IGVm = abs(alpha_0)-abs(teta_IGV_mid)/2-i_IGV_od;
gamma_IGVm =mod_gamma_IGVm;

mod_gamma_IGVt = abs(alpha_0)-abs(teta_IGV_tip)/2-i_IGV_od;
gamma_IGVt =mod_gamma_IGVt;

%-------------------------------------------------------------------------
%                              NACA A4K6

x_perc=[0 1.25 2.5 5 10 15 20 30 40 50 60 70 80 90 95 100];% chord percenage x/c
y_t_perc=[0 0.771 1.057 1.462 2.01 2.386 2.656 2.954 2.971 2.723 2.301 1.87 1.438 1.007 0.791 0];% half thickness t/c
y_c_perc=[0 0.792 1.357 2.248 3.531 4.42 5.04 5.71 5.82 5.516 4.891 4.011 2.922 1.642 0.912 0];% camber line y/c
dy_dx_perc=[Inf 0.5034 0.41 0.3131 0.2110 0.1483 0.1023 0.0359 0.0116 0.0478 0.0761 0.099 0.1184 0.1387 0.155 -Inf];% derivative

%------------------------------------------------------------------------

% ROTOR

c_IGV=chord;

x_IGV = x_perc*c_IGV;
y_t_IGV = c_IGV*y_t_perc*0.8;

% HUB 

cl_IGVh=teta_IGV_hub/25;

% equations

y_IGVh = c_IGV*y_c_perc*cl_IGVh;
dy_dx_IGVh = dy_dx_perc*cl_IGVh;

eps_IGVh = atand(dy_dx_IGVh);

xss_IGVh = x_IGV-y_t_IGV.*sind(eps_IGVh);
yss_IGVh = y_IGVh+y_t_IGV.*cosd(eps_IGVh);

xps_IGVh = x_IGV+y_t_IGV.*sind(eps_IGVh);
yps_IGVh = y_IGVh-y_t_IGV.*cosd(eps_IGVh);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVh = -(xss_IGVh*sind(gamma_IGVh)+yss_IGVh*cosd(gamma_IGVh));
AX_coord_ss_IGVh = xss_IGVh*cosd(gamma_IGVh)-yss_IGVh*sind(gamma_IGVh);

T_coord_ps_IGVh = -(xps_IGVh*sind(gamma_IGVh)+yps_IGVh*cosd(gamma_IGVh));
AX_coord_ps_IGVh = xps_IGVh*cosd(gamma_IGVh)-yps_IGVh*sind(gamma_IGVh);

% MID 

cl_IGVm=teta_IGV_mid/25;

% equations

y_IGVm = c_IGV*y_c_perc*cl_IGVm;
dy_dx_IGVm = dy_dx_perc*cl_IGVm;

eps_IGVm = atand(dy_dx_IGVm);

xss_IGVm = x_IGV-y_t_IGV.*sind(eps_IGVm);
yss_IGVm = y_IGVm+y_t_IGV.*cosd(eps_IGVm);

xps_IGVm = x_IGV+y_t_IGV.*sind(eps_IGVm);
yps_IGVm = y_IGVm-y_t_IGV.*cosd(eps_IGVm);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVm = -(xss_IGVm*sind(gamma_IGVm)+yss_IGVm*cosd(gamma_IGVm));
AX_coord_ss_IGVm = xss_IGVm*cosd(gamma_IGVm)-yss_IGVm*sind(gamma_IGVm);

T_coord_ps_IGVm = -(xps_IGVm*sind(gamma_IGVm)+yps_IGVm*cosd(gamma_IGVm));
AX_coord_ps_IGVm = xps_IGVm*cosd(gamma_IGVm)-yps_IGVm*sind(gamma_IGVm);

% TIP 

cl_IGVt=teta_IGV_tip/25;

% equations

y_IGVt = c_IGV*y_c_perc*cl_IGVt;
dy_dx_IGVt = dy_dx_perc*cl_IGVt;

eps_IGVt = atand(dy_dx_IGVt);

xss_IGVt = x_IGV-y_t_IGV.*sind(eps_IGVt);
yss_IGVt = y_IGVt+y_t_IGV.*cosd(eps_IGVt);

xps_IGVt = x_IGV+y_t_IGV.*sind(eps_IGVt);
yps_IGVt = y_IGVt-y_t_IGV.*cosd(eps_IGVt);

% Passage from x-y coordinates to T-ax

T_coord_ss_IGVt = -(xss_IGVt*sind(gamma_IGVt)+yss_IGVt*cosd(gamma_IGVt));
AX_coord_ss_IGVt = xss_IGVt*cosd(gamma_IGVt)-yss_IGVt*sind(gamma_IGVt);

T_coord_ps_IGVt = -(xps_IGVt*sind(gamma_IGVt)+yps_IGVt*cosd(gamma_IGVt));
AX_coord_ps_IGVt = xps_IGVt*cosd(gamma_IGVt)-yps_IGVt*sind(gamma_IGVt); 

% plot(xss_IGVh,yss_IGVh,xps_IGVh,yps_IGVh)

% plot(-T_coord_ss_IGVh,-AX_coord_ss_IGVh,'green',-T_coord_ps_IGVh,-AX_coord_ps_IGVh,'green')
% hold on

plot(-T_coord_ss_IGVm,-AX_coord_ss_IGVm,'red',-T_coord_ps_IGVm,-AX_coord_ps_IGVm,'red')
% hold on

% plot(-T_coord_ss_IGVt,-AX_coord_ss_IGVt,'black',-T_coord_ps_IGVt,-AX_coord_ps_IGVt,'black')

% title('Shape of rotor blade NACA A4K6 series in x-y plane in off design')
hold off
title('Blade profile at midspan in design/off-design')

% axis([-10 12 -2 12])

% All these last commented lines allow to draw the blade profiles for IGV
% at hub, mid and tip, as can be noticed only the mid section is not
% commented in order to plot the profiles at midspan in design and off-design
% condition to see the difference in blade pitch angle.

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%                          VELOCITY TRIANGLES
%------------------------------------------------------------------------
%------------------------------------------------------------------------

drawArrowin = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'black');  
drawArrowout = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'red');
drawArrowv3 = @(x,y) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,'green');

%                                ROTOR
%------------------------------------------------------------------------
% HUB

y_rh_w1 = [0 w_1ax];
y_rh_v1 = [0 v_1ax];
y_rh_u1 = [v_1ax v_1ax];
x_rh_w1 = [0 w_1th];
x_rh_v1 = [0 v_1th];
x_rh_u1 = [w_1th v_1th];

y_rh_w2 = [0 w_2ax];
y_rh_v2 = [0 v_2ax(end)];
y_rh_u2 = [v_2ax(end) v_2ax(end)];
x_rh_w2 = [0 w_2th];
x_rh_v2 = [0 v_2th];
x_rh_u2 = [w_2th v_2th];

y_rh_v3 = [0 v_3ax];
x_rh_v3 = [0 v_3th];

figure;
subplot(3,1,3);
drawArrowin(x_rh_w1,y_rh_w1); hold on
drawArrowin(x_rh_v1,y_rh_v1); hold on
drawArrowin(x_rh_u1,y_rh_u1); hold on
drawArrowout(x_rh_w2,y_rh_w2); hold on
drawArrowout(x_rh_v2,y_rh_v2); hold on
drawArrowout(x_rh_u2,y_rh_u2); hold on
drawArrowv3(x_rh_v3,y_rh_v3); hold off
axis ij
title('Velocity triangle at rotor hub')
xlabel('Tangential direction')
ylabel('Axial direction')

% MID

y_rm_w1 = [0 w_1ax];
y_rm_v1 = [0 v_1ax];
y_rm_u1 = [v_1ax v_1ax];
x_rm_w1 = [0 w_1tm];
x_rm_v1 = [0 v_1tm];
x_rm_u1 = [w_1tm v_1tm];

y_rm_w2 = [0 w_2ax];
y_rm_v2 = [0 v_2ax(end)];
y_rm_u2 = [v_2ax(end) v_2ax(end)];
x_rm_w2 = [0 w_2tm];
x_rm_v2 = [0 v_2tm];
x_rm_u2 = [w_2tm v_2tm];

y_rm_v3 = [0 v_3ax];
x_rm_v3 = [0 v_3tm];

subplot(3,1,2);
drawArrowin(x_rm_w1,y_rm_w1); hold on
drawArrowin(x_rm_v1,y_rm_v1); hold on
drawArrowin(x_rm_u1,y_rm_u1); hold on
drawArrowout(x_rm_w2,y_rm_w2); hold on
drawArrowout(x_rm_v2,y_rm_v2); hold on
drawArrowout(x_rm_u2,y_rm_u2); hold on
drawArrowv3(x_rm_v3,y_rm_v3); hold off
axis ij
title('Velocity triangle at rotor mid')
xlabel('Tangential direction')
ylabel('Axial direction')


% TIP

y_rt_w1 = [0 w_1ax];
y_rt_v1 = [0 v_1ax];
y_rt_u1 = [v_1ax v_1ax];
x_rt_w1 = [0 w_1tt];
x_rt_v1 = [0 v_1tt];
x_rt_u1 = [w_1tt v_1tt];

y_rt_w2 = [0 w_2ax];
y_rt_v2 = [0 v_2ax(end)];
y_rt_u2 = [v_2ax(end) v_2ax(end)];
x_rt_w2 = [0 w_2tt];
x_rt_v2 = [0 v_2tt];
x_rt_u2 = [w_2tt v_2tt];

y_rt_v3 = [0 v_3ax];
x_rt_v3 = [0 v_3tt];

subplot(3,1,1);
drawArrowin(x_rt_w1,y_rt_w1); hold on
drawArrowin(x_rt_v1,y_rt_v1); hold on
drawArrowin(x_rt_u1,y_rt_u1); hold on
drawArrowout(x_rt_w2,y_rt_w2); hold on
drawArrowout(x_rt_v2,y_rt_v2); hold on
drawArrowout(x_rt_u2,y_rt_u2); hold on
drawArrowv3(x_rt_v3,y_rt_v3); hold off
axis ij
title('Velocity triangle at rotor tip')
xlabel('Tangential direction')
ylabel('Axial direction')

%------------------------------------------------------------------------
% Velocity triangles in each section

figure;
drawArrowin(x_rh_w1,y_rh_w1); hold on
drawArrowin(x_rh_v1,y_rh_v1); hold on
drawArrowin(x_rh_u1,y_rh_u1); hold on
drawArrowout(x_rh_w2,y_rh_w2); hold on
drawArrowout(x_rh_v2,y_rh_v2); hold on
drawArrowout(x_rh_u2,y_rh_u2); hold on
drawArrowv3(x_rh_v3,y_rh_v3); hold off
axis ij
title('Velocity triangle at rotor hub')
xlabel('Tangential direction')
ylabel('Axial direction')

figure;
drawArrowin(x_rm_w1,y_rm_w1); hold on
drawArrowin(x_rm_v1,y_rm_v1); hold on
drawArrowin(x_rm_u1,y_rm_u1); hold on
drawArrowout(x_rm_w2,y_rm_w2); hold on
drawArrowout(x_rm_v2,y_rm_v2); hold on
drawArrowout(x_rm_u2,y_rm_u2); hold on
drawArrowv3(x_rm_v3,y_rm_v3); hold off
axis ij
title('Velocity triangle at rotor mid')
xlabel('Tangential direction')
ylabel('Axial direction')

figure;
drawArrowin(x_rt_w1,y_rt_w1); hold on
drawArrowin(x_rt_v1,y_rt_v1); hold on
drawArrowin(x_rt_u1,y_rt_u1); hold on
drawArrowout(x_rt_w2,y_rt_w2); hold on
drawArrowout(x_rt_v2,y_rt_v2); hold on
drawArrowout(x_rt_u2,y_rt_u2); hold on
drawArrowv3(x_rt_v3,y_rt_v3); hold off
axis ij
title('Velocity triangle at rotor tip')
xlabel('Tangential direction')
ylabel('Axial direction')







