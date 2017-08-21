%script to creat movie for protein concentration

function M = GFP_movie(Nx, fname)

RD;

LN = length(L);  %numver of saved results

Lmax = max(L);  % maximum nondimensional biofilm thickness

Wmax = max(max(W));

%Nx = 100;

for i = 1:LN
    
    hold off;
    
    dx = L(i)/Nx;
    
    x = (dx:dx:L(i)) - .5*dx; % coordinate of offset grid
    
    hw = area(x,W(i,:));
    
    set(hw,'FaceColor',[0 max(W(i,:))/Wmax 0]);
    
    hold on;
    
    axis([0 1.1*Lmax 0 1.1*Wmax]);
    
    yE = 0:.1*Wmax:1.1*Wmax;
    xE = L(i)*ones(size(yE)); % x-coordinate of biofilm-bulk interface
    
    plot(xE,yE,'r-','linewidth',2);
    
    xlabel('x (unit = 139 \mu m)','fontsize',14)
    ylabel('W_i','fontsize',14);
    
    tL = sprintf('t = %d hours',i-1);
    text(0.5, 0.9, tL,'fontsize',14);
    
    %pause(.3);
    
    M(i) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);

end

movie2avi(M,fname,'FPS',6,'compression','None');
    
    
    
    
    
    
    

