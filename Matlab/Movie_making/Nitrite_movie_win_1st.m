%script to creat movie for protein concentration

function M = Nitrite_movie_win_1st(L, Nx, fname)

RD;

Wmax = max(max(W));
SO2max = max(SO2);
SNO2max = max(SNO2);

Allmax = max([Wmax SO2max SNO2max]);

%Nx = 100;
figure(1);



hold on;

dx = L/Nx;

x = 0:dx:L;


[Wi, Wj] = size(W);

for i = 1:Wi
    hold off;
        
    hw = area(x,W(i,:));
   
    set(hw,'FaceColor',[0 max(W(i,:))/Wmax 0]);
    hold on;
    plot(x,SO2,'r-','linewidth',2);
    axis([0 L 0 1.1*Allmax]);
  
    plot(x,SNO2,'-','linewidth',2);
    legend('mRNA','Oxygen','Nitrite','Location','Northwest');
    

        
%     axis([0 1.1*Lmax 0 1.1*Wmax]);
    
%     yE = 0:.1*Wmax:1.1*Wmax;
%     xE = L(i)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
%     
%     plot(xE,yE,'b-','linewidth',2);
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('W_i','fontsize',14);
    
    tL = sprintf('t = %d hours',i);
    text(0.5, 0.6, tL,'fontsize',14);
    
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
    pause;
end

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    

