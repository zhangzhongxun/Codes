%plot script

dNx = 5;
[Wi,Wj] = size(QS);
dx = 1/100;
x = .5*dx:dx:1-.5*dx;
y = 0:72;
[X,Y] = meshgrid(x,y);
surf(X,Y,QS(:,dNx:dNx:Wj));
%surf(X,Y,WL);
%surf(X,Y,QS);
%surf(X,Y,S);
xlabel('\zeta','fontsize',14)
ylabel('t (hours)','fontsize',14)
zlabel('QS','fontsize',14)