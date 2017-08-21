%Script to test codes for Acid Stress inducible gene

Set_ParaA;

global dx Nx dt t0 PD h0

x = 0:dx:1;

p = 2*(x-1).^2;

r = 1 - p;

r = 12.85*r(1:Nx);

% [d, e, f, rhs] = Build_SA(r);
%     
% S = Tridiag_Solver(d,e,f, rhs);

% mu = (0.8/3600)*(1 - p/6);
% 
% f_asr = ones(size(mu));
% pos = find(p < 1);
% f_asr(pos) = p(pos).^3;
% r = t0*f_asr.*mu;

% wn = 1 - x.^2;

% L = 1.2;
% 
% [L_dot, u] = Velo_growthA(L, mu);
% 
% 
% 
% 
% 
% [d, e, f, rhs] = Build_WA(r, mu, wn, L_dot, L, u);
%     
% WA = Tridiag_Solver(d,e,f, rhs);

L = 1.4;

S = Solve_SAPD(r,L);

% fid1 = fopen('./data/W','w'); 
% 
% fprintf(fid1, '%12.10e ', WA);
% fprintf(fid1, '\n');
% fclose(fid1);

