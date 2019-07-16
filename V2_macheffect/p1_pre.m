  V_0T_m = 0;
  V_0T_t = 0; 
  V_0T_h = 0;
  
  % Hub diameter
  D_h = D_t * lambda;
  b =  (D_t - D_h)/2 ;
  D_m = (D_h + D_t)/2;
  U_m = omega * D_m / 2;
    
  c_IGV_t = c_IGV;
  c_IGV_h = c_IGV;