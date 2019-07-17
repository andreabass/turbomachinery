function Y = y_AM_od(alpha_in_geo,alpha_in,alpha_out,sigma,inc,rho_out,b,c,V_out,mu)


beta_1prime_od = 90-alpha_in_geo;
alpha_2prime_od = 90-alpha_out;

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

Xi_AM = 1/sigma-0.75;

if 1/sigma <= 0.8
    Dis = -38*Xi_AM-53.5*Xi_AM^2-29*Xi_AM^3;
else
    Dis = 2.0374-(1/sigma-0.8)*(69.58-(alpha_2prime_od/14.48)^3.1);
end

is_AM = i_sr+Dis;

if inc/is_AM < -3
    Kinc = -1.39214-1.90738*(inc/is_AM);
elseif inc/is_AM >=-3 && inc/is_AM <0
    Kinc = 1+0.52*(abs(inc/is_AM))^1.7;
elseif inc/is_AM >=0 && inc/is_AM <1.7
    Kinc = 1+(inc/is_AM)^(2.3+0.5*inc/is_AM);
else
    Kinc = 6.23+9.8577*(inc/is_AM-1.7);
end
             
C = 0.08*((alpha_2prime_od/30)^2-1);
n_AM = 1+alpha_2prime_od/30;

if alpha_2prime_od <= 27
    A = 0.025+(27-alpha_2prime_od)/530;
else
    A = 0.025+(27-alpha_2prime_od)/3850;
end

if alpha_2prime_od <= 30
    s_over_c_min = 0.46+alpha_2prime_od/77; % ratio s/c which minimize the profile losses
    B = 0.1583-alpha_2prime_od/1640;
    X_AM = 1/sigma-s_over_c_min;
    Y_inc = A+B*X_AM^2+C*X_AM^3; 
else 
    s_over_c_min = 0.614+alpha_2prime_od/130; % ratio s/c which minimize the profile losses
    B = 0;
    X_AM = 1/sigma-s_over_c_min;
    Y_inc = A+B*(abs(X_AM))^n_AM;
end 

    Y_inc = Y_inc*Kinc; 
    
    Re_out = rho_out * V_out * c / mu;
    alpha_av_t_01 = atand((tand(alpha_in)+tand(alpha_out))/2);
    cL = 2 * s_over_c_min * (abs(tand(alpha_out)-tand(alpha_in)))*cosd(alpha_av_t_01);

        Y_p_1_Re = Y_inc * (2e5 / Re_out)^0.2;
        Y_1_sec = c / b *(0.0334 * cosd(alpha_out)/cosd(alpha_in)) * (cL / s_over_c_min)^2 * ((cosd(alpha_out))^2) / (cosd(alpha_av_t_01))^3;
    
    Y = Y_p_1_Re + Y_1_sec;

end

