%script to creat movie for protein concentration

function M = S_GFP_movie_win(Nx, fname)

RD;

LN = length(L);  %numver of saved results

Lmax = max(L);  % maximum nondimensional biofilm thickness

Wmax = max(max(W));
Smax = max(max(S));

%Nx = 100;
figure(1);

for i = 1:LN
    
    
    
    dx = L(i)/Nx;
    
    x = (dx:dx:L(i)) - .5*dx; % coordinate of offset grid
    
    subplot(1,2,1);
    hold off;
    
    hs = area(x,S(i,:));
    set(hs,'BaseValue',min(S(i,:)));
    set(hs,'FaceColor',[max(S(i,:))/Smax 0 0]);
  
    hold on;
    axis([0 1.1*Lmax 0 1.1*Smax]);
    
    yE = 0:.1*Smax:1.1*Smax;
    xE = L(i)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'b-','linewidth',2);
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('S','fontsize',14);
    
    tL = sprintf('t = %d hours',i-1);
    text(0.5, 0.9, tL,'fontsize',14);
    
    subplot(1,2,2);
    hold off;
    
    hw = area(x,W(i,:));
    
    set(hw,'FaceColor',[0 max(W(i,:))/Wmax 0]);
    
    hold on;
    
    axis([0 1.1*Lmax 0 1.1*Wmax]);
    
    yE = 0:.1*Wmax:1.1*Wmax;
    xE = L(i)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'b-','linewidth',2);
    
    xlabel('x (unit = 100 \mu m)','fontsize',14)
    ylabel('W_i','fontsize',14);
    
    tL = sprintf('t = %d hours',i-1);
    text(0.5, 0.9, tL,'fontsize',14);
    
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
end

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    
