  V_0T_m = 0;
  V_0T_t = 0; 
  V_0T_h = 0;
    
  phi_1_m = 0.5;
          
  b =  (D_t - D_h(end))/2 ;
  D_m = (D_h(end) + D_t)/2;
  U_m = omega * D_m / 2;
    
  V_1A_m     = phi_1_m * omega * D_m / 2;
          V_1A       = V_1A_m; % (definition)
    
  V_1A_t = V_1A_m;
  V_1A_h = V_1A_m;
    
  c_IGV_t = c_IGV;
  c_IGV_h = c_IGV;