%plot script
dNx = 5;
[Wi,Wj] = size(WL);

dx = 1/100;
x = .5*dx:dx:1-.5*dx;
y = 0:72;
[X,Y] = meshgrid(x,y);
surf(X,Y,WL(:,dNx:dNx:Wj));
%surf(X,Y,WL);
%surf(X,Y,QS);
%surf(X,Y,S);
xlabel('\zeta','fontsize',20)
ylabel('t (hours)','fontsize',20)
zlabel('lsaI','fontsize',20)