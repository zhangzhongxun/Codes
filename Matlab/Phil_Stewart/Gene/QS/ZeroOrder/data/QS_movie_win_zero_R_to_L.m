%script to creat movie for protein concentration

function M = QS_movie_win_zero_R_to_L(Nx, fname)

RD;

%  QS = QS(1:49,5:5:500);
%  S = S(1:49,5:5:500);
 
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

subplot(1,1,1);
tmpmax = max(Smax,QSmax);

[Wi, Wj] = size(WL);

for i = 1:LN
    
    dx = L(i)/Nx;
    
    x = (dx:dx:L(i)) - .5*dx;
    
    x = Lmax - x;
    
    hold off;
    
    hw = area(x,QS(i,:));
   
    set(hw,'FaceColor',[0 max(QS(i,:))/QSmax 0]);
    hold on;
    plot(x,S(i,:),'r-','linewidth',2);
    axis([0 Lmax-0.5*dx 0 1.1*tmpmax]);
    
    yE = 0:.1*tmpmax:1.1*tmpmax;
    xE = (Lmax-L(i)+.5*dx)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'b-','linewidth',2);
    legend('QS','Oxygen');
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('QS','fontsize',14);
    
    tL = sprintf('t = %d hours',i);
    text(.5*Lmax, 0.8*tmpmax, tL,'fontsize',14);
    
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
end

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    

