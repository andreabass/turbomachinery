function velt(V,W,U,varargin)

if isempty(varargin)
    style = 'k';
else style = varargin{1};
end

if U==0
    
hold on
quiver(0,0,W,V,0,style,'LineWidth',2,'MaxHeadSize',0.15);
    
else
    
beta = acosd( (W^2+U^2-V^2) / (2*W*U) )  - 90;
Wa   = W * cosd(beta);
Wt   = W * sind(beta);
Va   = Wa;
Vt   = Wt + U;

hold on
quiver(0,0,Vt,Va,0,style,'LineWidth',2,'MaxHeadSize',0.15);
quiver(0,0,Wt,Wa,0,style,'LineWidth',2,'MaxHeadSize',0.15);
quiver(Wt,Wa,U,0,0,style,'LineWidth',2,'MaxHeadSize',0.15);

end

YYend = 300;
YYstart = 0;
XXend = 300;
XXstart = -300;

axis([XXstart XXend YYstart YYend])
pbaspect([1 (YYend-YYstart)/(XXend-XXstart) 1])
axis ij
xlabel('Tangential direction')
ylabel('Axial direction')

end

