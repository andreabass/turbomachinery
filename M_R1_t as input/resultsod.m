p_T3_av     = sum(p_T3_ave.*dm)/m;
beta_TT     = p_T3_av/p_T0;
T_T3_av     = sum(T_T3_ave.*dm)/m;
deltaHis_TT = cp * T_T0 * (beta_TT ^ ((gamma - 1)/gamma)-1); % [J/kg]
eta_TT      =  deltaHis_TT / ( cp*(T_T3_av - T_T0) );

Yigv    = (p_T0 - p_T1) ./ (p_T1 - p_1);
Yigv_av = sum(Yigv)/length(Yigv);

Yrot    = (p_TR1 - p_TR2) ./ (p_TR1 - p_1);
Yrot_av = sum(Yrot)/length(Yrot);

Ystat    = (p_T2 - p_T3) ./ (p_T2 - p_2);
Ystat_av = sum(Ystat)/length(Ystat);

[Yigv_av Yrot_av Ystat_av] * 100

