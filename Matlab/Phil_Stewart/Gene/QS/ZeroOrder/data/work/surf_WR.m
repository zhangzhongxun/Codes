%plot script
dNx = 5;
[Wi,Wj] = size(WR);

dx = 1/100;
x = .5*dx:dx:1-.5*dx;
y = 0:72;
[X,Y] = meshgrid(x,y);
surf(X,Y,WR(:,dNx:dNx:Wj));
%surf(X,Y,WL);
%surf(X,Y,QS);
%surf(X,Y,S);
xlabel('\zeta','fontsize',20)
ylabel('t (hours)','fontsize',20)
zlabel('rsaL','fontsize',20)