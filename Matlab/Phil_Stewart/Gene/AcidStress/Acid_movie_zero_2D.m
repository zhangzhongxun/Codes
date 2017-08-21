%script to creat movie for protein concentration

function M = Acid_movie_zero_2D(Nx, fname, flag)

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
close(1);
figure(1);

tmpmax = max(Smax,Pmax);

[Wi, Wj] = size(W);

for i = 2:8  % 2:LN
    
    dr = L(i)/Nx;
    
    r = 0:dr:L(i);
    
    theta = 0:pi/30:pi;
    
    [R,T] = meshgrid(r,theta);
    
    X = R.*cos(T);
    Y = R.*sin(T);
    
    [Si, tmp] = meshgrid(S(i,:),theta);
    [Pi, tmp] = meshgrid(P(i,:),theta);
    [Wi, tmp] = meshgrid(W(i,:),theta);
    
    Smi = max(S(i,:));
    Pmi = max(P(i,:));
    Wmi = max(W(i,:));
    
    if flag == 1 % plot glucose
        
        contour(X,Y,R);
        
        %plot glucose
        % subplot(1,3,1);
%         hold off;
%         cmr = zeros(64,3);
%         cmr(:,1) = 0:(Smi/Smax)/63:(Smi/Smax)';
%         
%         contourf(X,Y,Si,40);
%         colormap(cmr);
%         colorbar;
%         %freezeColors;
%         %cbfreeze(colorbar);
%         shading flat;
%         hold on;
%         
%         axis([-1.1*Lmax 1.1*Lmax 0 1.1*Lmax]);
%         axis equal;
%         
%         xlabel('x (unit = 100 \mu m)','fontsize',14)
%         ylabel('Normalized Glucose Concentration','fontsize',14);
%         
%         %cbfreeze;
%         
%         tL = sprintf('t = %d hours',i);
%         text(-.2*Lmax, 1.3*Lmax, tL,'fontsize',14);
        
    elseif flag == 2 % plot lactate
        % plot lactate
        %subplot(1,3,2);
        hold off;
        cmb = zeros(64,3);
        cmb(:,3) = 0:(Pmi/Pmax)/63:(Pmi/Pmax)';
        
        contourf(X,Y,Pi,40);
        colormap(cmb);
        colorbar;
        %     freezeColors;
        %     cbfreeze(colorbar);
        shading flat;
        hold on;
        
        axis([-1.1*Lmax 1.1*Lmax 0 1.1*Lmax]);
        axis equal;
        
        xlabel('x (unit = 100 \mu m)','fontsize',14)
        ylabel('Normalized Lactate Concentration','fontsize',14);
        colorbar;
        %cbfreeze;
        
        tL = sprintf('t = %d hours',i);
        text(-.2*Lmax, 1.3*Lmax, tL,'fontsize',14);
        
    elseif flag == 3 % plot gene
        %plot acid stress gene
        %subplot(1,3,3);
        hold off;
        cmg = zeros(64,3);
        cmg(:,2) = 0:(Wmi/Wmax)/63:(Wmi/Wmax)';
        
        contourf(X,Y,Wi,40);
        colormap(cmg);
        colorbar;
        %     freezeColors;
        %     cbfreeze(colorbar);
        shading flat;
        hold on;
        
        axis([-1.1*Lmax 1.1*Lmax 0 1.1*Lmax]);
        axis equal;
        
        xlabel('x (unit = 100 \mu m)','fontsize',14)
        ylabel('Normalized Acid Stress Gene Concentration','fontsize',14);
        
        %cbfreeze;
        
        tL = sprintf('t = %d hours',i);
        text(-.2*Lmax, 1.3*Lmax, tL,'fontsize',14);
        
    end
    
    M(i-1) = getframe(gcf);
    fprintf(1,'frame %d done\n', i);
    
end

fprintf(1,'It is ok until here\n');
movie(M);

movie2avi(M,fname,'FPS',6,'compression','Cinepak');
    
    
    
    
    
    
    

