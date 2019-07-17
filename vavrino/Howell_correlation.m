load('Howell_1');
load('Howell_2');

Dbeta_Psi_curve = spline(Howell_1(:,1),Howell_1(:,2));
Psi_curve = spline(Howell_2(:,1),Howell_2(:,2));
