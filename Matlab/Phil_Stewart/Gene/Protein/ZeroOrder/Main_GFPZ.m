Set_ParaZ;   % Set the parameter values for simulation

global  Mu_max Rho Dw Ds Sigma LSc Sc % global parameters, most of them are dimensionless

global Nx dx dt FT FS_interval flag_w flag_kg% control parameters

% Open files to save computed data, it is saved every FS_interval time
% steps


fid1 = fopen('./data/W','w');  % file for W :GFP

fid2 = fopen('./data/S','w');      % file for S: substrate

fid3 = fopen('./data/u','w');        % file for u: velocity

fid4 = fopen('./data/L','w');        % file for L: biofilm thickness

[L, w_n, s, u] = Initial_GFPZ(); 

% Save the initial condition

fprintf(fid1, '%12.10e ', w_n);
fprintf(fid1, '\n');

fprintf(fid2, '%12.10e ', s);
fprintf(fid2, '\n');

fprintf(fid3, '%12.10e ', u);
fprintf(fid3, '\n');

fprintf(fid4, '%12.10e ', L);
fprintf(fid4, '\n');

for t_number = 1: FT/dt
    
    t_current = t_number*dt;
    
    if mod(t_number,FS_interval) == 0
        
        fprintf('\n Time = %6.4f, Time steps = %d \n', t_current, t_number);
        
    end
    
%     [d, e, f, rhs] = Build_SZ(L, 1);
%     
%     s = Tridiag_Solver(d,e,f, rhs);

    s = Solve_SZ(L, 1);
    
    [L_dot, u] = Velo_growthZ(L, s);
    
    
    fi = 1;  % turn on inducer after 48 hours, simulation starts now
    
    [d, e, f, rhs] = Build_WZ(L, L_dot, s, w_n, fi, u);
    
    w_n = Tridiag_Solver(d,e,f, rhs);
    
    L = L + L_dot*dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Save intermediate results if t_number is a multiple of FS_interval     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(t_number,FS_interval) == 0
        fprintf(fid1, '%12.10e ', w_n);
        fprintf(fid1, '\n');
        
        fprintf(fid2, '%12.10e ', s);
        fprintf(fid2, '\n');
        
        fprintf(fid3, '%12.10e ', u);
        fprintf(fid3, '\n');
        
        fprintf(fid4, '%12.10e ', L);
        fprintf(fid4, '\n');
        
    end
    
end
    
    
    
    