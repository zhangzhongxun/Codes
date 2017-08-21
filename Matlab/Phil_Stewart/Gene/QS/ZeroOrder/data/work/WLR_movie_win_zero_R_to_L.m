%script to creat movie for protein concentration

function M = WLR_movie_win_zero_R_to_L(Nx,dNx, fname)

RD;

 WL = WL(1:49,5:5:500);
 WR = WR(1:49,5:5:500);
 S = S(1:49,5:5:500);
 
% L = L(1:49);
 
LN = length(L);  %numver of saved results

Lmax = max(L); 

QSmax = max(max(QS));
WLmax = max(max(WL));
WRmax = max(max(WR));
Smax = max(max(S));

Allmax = max([WLmax WRmax Smax QSmax]);

%Nx = 100;
figure(1);


SLmax = max(Smax,WLmax);
SRmax = max(Smax,WRmax);

[Wi, Wj] = size(WL);

for i = 1:LN
    
    dx = L(i)/Nx;
    
    x = (dx:dx:L(i)) - .5*dx;
    
    x = Lmax - x;
    
    subplot(1,2,1);
    hold off;
    
    hw = area(x,WL(i,:));
   
    set(hw,'FaceColor',[0 max(WL(i,:))/WLmax 0]);
    hold on;
    plot(x,S(i,:),'r-','linewidth',2);
    axis([0 Lmax-0.5*dx 0 1.1*SLmax]);
    
    yE = 0:.1*SLmax:1.1*SLmax;
    xE = (Lmax-L(i)+.5*dx)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'b-','linewidth',2);
    legend('lasI','Oxygen');
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('lasI','fontsize',14);
    
    tL = sprintf('t = %d hours',i);
    text(.5*Lmax, 0.8*SLmax, tL,'fontsize',14);
    
    subplot(1,2,2);
    hold off;
    
    hw = area(x,WR(i,:));
   
    set(hw,'FaceColor',[0 0 max(WR(i,:))/WRmax]);
    hold on;
    plot(x,S(i,:),'r-','linewidth',2);
    axis([0 Lmax-0.5*dx 0 1.1*SRmax]);
    
    yE = 0:.1*SRmax:1.1*SRmax;
    xE = (Lmax-L(i)+.5*dx)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'b-','linewidth',2);
    legend('rsaL','Oxygen');
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('rsaL','fontsize',14);
    
    tL = sprintf('t = %d hours',i);
    text(.5*Lmax, 0.8*SRmax, tL,'fontsize',14);
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
end

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    

