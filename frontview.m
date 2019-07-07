figure

xouter = linspace(-D_t/2,D_t/2,1000);
ypouter = sqrt((D_t/2)^2-xouter.^2);
ymouter = -sqrt((D_t/2)^2-xouter.^2);
xinner = linspace(-D_h/2,D_h/2,1000);
ypinner = sqrt((D_h/2)^2-xinner.^2);
yminner = -sqrt((D_h/2)^2-xinner.^2);
plot(xouter,ypouter,'k',xouter,ymouter,'k','LineWidth',3)
hold on
axis([ -0.7 0.7 -0.7 0.7])
pbaspect([1 1 1])

plot(xinner,ypinner,'k',xinner,yminner,'k','LineWidth',3)

angle = 2*pi/N_R;
angles = 0;
for i = 2:N_R
angles(i) = angles(i-1)+angle;
end

coordxouter = D_t/2*cos(angles);
coordyouter = D_t/2*sin(angles);

coordxinner = D_h/2*cos(angles);
coordyinner = D_h/2*sin(angles);

 for i = 1:N_R
     plot([coordxinner(i) coordxouter(i)], [coordyinner(i) coordyouter(i)], 'k', 'LineWidth',3)
 end
 
 title('FRONT VIEW (ROTOR)')