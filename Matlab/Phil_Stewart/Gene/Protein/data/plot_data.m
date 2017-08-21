
plot(xm,WND(51,:),'g-o')
hold on;
plot(xm,WND(81,:),'r-.h')
plot(xm,WND(101,:),'-->')
xlabel('\zeta','fontsize',14)
ylabel('W_i','fontsize',14)
legend('t = 50 h','t = 80 h','t = 100 h')
subplot(1,2,2)
plot(xm,WD(51,:),'g-o')
hold on;
plot(xm,WD(81,:),'r-.h')
plot(xm,WD(101,:),'-->')
xlabel('\zeta','fontsize',14)
ylabel('W_i','fontsize',14)
legend('t = 50 h','t = 80 h','t = 100 h')
figure(3)
