clear all;

Set_ParaA;   % Set the parameter values for simulation

%global  Mu_max  Dq Ds Sigma LSc Sc % global parameters, most of them are dimensionless

global  dt FT FS_interval G_s G_p Rho t0 % control parameters

% Open files to save computed data, it is saved every FS_interval time
% steps

fid1 = fopen('./data/S','w');  % file for S: glucose

fid2 = fopen('./data/P','w');  % file for P: lactate

fid3 = fopen('./data/W','w');  % file for W: gene

fid4 = fopen('./data/u','w');        % file for u: velocity

fid5 = fopen('./data/L','w');        % file for L: biofilm thickness

fid6 = fopen('./data/mu','w');        % file for mug: growth-rate

[L, Wn, P, S, u] = Initial_Acid();
% Save the initial condition
mug = mu_growth(P, S);

fprintf(fid1, '%12.10e ', S);
fprintf(fid1, '\n');

fprintf(fid2, '%12.10e ', P);
fprintf(fid2, '\n');

fprintf(fid3, '%12.10e ', Wn);
fprintf(fid3, '\n');

fprintf(fid4, '%12.10e ', u);
fprintf(fid4, '\n');

fprintf(fid5, '%12.10e ', L);
fprintf(fid5, '\n');

fprintf(fid6, '%12.10e ', mug);
fprintf(fid6, '\n');

for t_number = 1: FT/dt
    
    t_current = t_number*dt;
    
    if mod(t_number,FS_interval) == 0
        
        fprintf('\n Time = %6.4f, Time steps = %d \n', t_current, t_number);
        
    end
    
    % calculate the current biofilm growth rate mu
    
    %r_S = (L^2)*G_s*mug;
    
    S = Solve_SAPD(P, L);  % calculate the glucose concentration
    
    mug = mu_growth(P, S);
    
    r_P = (L^2)*G_p*mug;
    
    P = Solve_PA(r_P);  % calculate the lactate concentration
    
    [L_dot, u] = Velo_growthA(L, mug);
    
    fasr = f_asr(P);
    
    r_W = Rho*t0*fasr.*mug;
    
    [d, e, f, rhs] = Build_WA(r_W, mug, Wn, L_dot, L, u);
    
    Wn = Tridiag_Solver(d,e,f, rhs);
    
    L = L + L_dot*dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Save intermediate results if t_number is a multiple of FS_interval     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(t_number,FS_interval) == 0
      
        fprintf(fid1, '%12.10e ', S);
        fprintf(fid1, '\n');
        
        fprintf(fid2, '%12.10e ', P);
        fprintf(fid2, '\n');
        
        fprintf(fid3, '%12.10e ', Wn);
        fprintf(fid3, '\n');
        
        fprintf(fid4, '%12.10e ', u);
        fprintf(fid4, '\n');
        
        fprintf(fid5, '%12.10e ', L);
        fprintf(fid5, '\n');
        
        fprintf(fid6, '%12.10e ', mug);
        fprintf(fid6, '\n');
        
    end
    
end
    
    fclose('all');
    
    