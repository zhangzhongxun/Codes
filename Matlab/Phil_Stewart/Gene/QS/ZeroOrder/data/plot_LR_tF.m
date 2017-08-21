%plot script

dNx = 5;
[Wi,Wj] = size(WL);

hold off;
dx = 1/100;
x = .5*dx:dx:1-.5*dx;
plot(x,WL(73,dNx:dNx:Wj),'-o');
hold on;
plot(x,WR(73,dNx:dNx:Wj),'g-x');
xlabel('\zeta','fontsize',14);
ylabel('Profile at t = 72 hours','fontsize',14);
legend('lsaI','rsaL');