Set_ParaQS;   % Set the parameter values for simulation

%global  Mu_max  Dq Ds Sigma LSc Sc % global parameters, most of them are dimensionless

global  dt FT FS_interval Q_BC  % control parameters

% Open files to save computed data, it is saved every FS_interval time
% steps

fid1 = fopen('./data/QS','w');  % file for QS: Quorum sensing molecues

fid2 = fopen('./data/WL','w');  % file for WL : lsaI

fid3 = fopen('./data/WR','w');  % file for WR : rsaL

fid4 = fopen('./data/S','w');      % file for S: substrate

fid5 = fopen('./data/u','w');        % file for u: velocity

fid6 = fopen('./data/L','w');        % file for L: biofilm thickness

[L, WL_n, WR_n, QS_n, S, u] = Initial_QS();
% Save the initial condition

fprintf(fid1, '%12.10e ', QS_n);
fprintf(fid1, '\n');

fprintf(fid2, '%12.10e ', WL_n);
fprintf(fid2, '\n');

fprintf(fid3, '%12.10e ', WR_n);
fprintf(fid3, '\n');

fprintf(fid4, '%12.10e ', S);
fprintf(fid4, '\n');

fprintf(fid5, '%12.10e ', u);
fprintf(fid5, '\n');

fprintf(fid6, '%12.10e ', L);
fprintf(fid6, '\n');


for t_number = 1: FT/dt
    
    t_current = t_number*dt;
    
    if mod(t_number,FS_interval) == 0
        
        fprintf('\n Time = %6.4f, Time steps = %d \n', t_current, t_number);
        
    end
    
    [S, x_S_0] = Solve_S(L, 1);
    
    [L_dot, u] = Velo_growth(L, S, x_S_0);
    
    [d, e, f, rhs] = Build_QS(L, L_dot, S, WL_n, u, Q_BC);
    
    QS_n = Tridiag_Solver(d,e,f, rhs);
    
    [d, e, f, rhs] = Build_WL(L, L_dot, S, WL_n, WR_n, QS_n, u);
    
    WL_n = Tridiag_Solver(d,e,f, rhs);
    
    [d, e, f, rhs] = Build_WR(L, L_dot, S, WR_n, QS_n, u);
    
    WR_n = Tridiag_Solver(d,e,f, rhs);
    
    L = L + L_dot*dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Save intermediate results if t_number is a multiple of FS_interval     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(t_number,FS_interval) == 0
        fprintf(fid1, '%12.10e ', QS_n);
        fprintf(fid1, '\n');
        
        fprintf(fid2, '%12.10e ', WL_n);
        fprintf(fid2, '\n');
        
        fprintf(fid3, '%12.10e ', WR_n);
        fprintf(fid3, '\n');
        
        fprintf(fid4, '%12.10e ', S);
        fprintf(fid4, '\n');
        
        fprintf(fid5, '%12.10e ', u);
        fprintf(fid5, '\n');
        
        fprintf(fid6, '%12.10e ', L);
        fprintf(fid6, '\n');
        
    end
    
end
    
    fclose('all');
    
    