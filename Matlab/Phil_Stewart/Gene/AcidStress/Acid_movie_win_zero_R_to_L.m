%script to creat movie for protein concentration

function M = Acid_movie_win_zero_R_to_L(Nx, fname)

RD;

%  QS = QS(1:49,5:5:500);
%  S = S(1:49,5:5:500);
 
% L = L(1:49);
 
LN = length(L);  %numver of saved results

Lmax = max(L); 

Smax = max(max(S));
Wmax = max(max(W));
Pmax = max(max(P));

Allmax = max([Wmax Smax Pmax]);

%Nx = 100;
figure(1);

subplot(1,1,1);
tmpmax = max(Smax,Pmax);

[Wi, Wj] = size(W);

for i = 1:LN
    
    dx = L(i)/Nx;
    
    x = 0:dx:L(i);
    
    x = Lmax - x;
    
    hold off;
    
    hw = area(x,W(i,:));
   
    set(hw,'FaceColor',[0 1 0]);
    hold on;
    plot(x,S(i,:),'r-','linewidth',2);
    axis([-.1*Lmax Lmax 0 1.1*Allmax]);
    plot(x,P(i,:),'-','linewidth',2);
   % axis([0 Lmax-0.5*dx 0 1.1*tmpmax]);
    
    
    
    yE = 0:.1*Allmax:1.1*Allmax;
    xE = (Lmax-L(i))*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'k-','linewidth',2);
   legend('Wi','Glucose','Lactate','Location','Northeast');
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('Acid Stress Gene','fontsize',14);
    
    tL = sprintf('t = %d hours',i);
    text(.5*Lmax, 1.0*Allmax, tL,'fontsize',14);
    
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
    pause;
end

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    

