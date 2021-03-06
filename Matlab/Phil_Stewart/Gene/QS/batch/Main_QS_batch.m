Set_ParaQS;   % Set the parameter values for simulation

%global  Mu_max  Dq Ds Sigma LSc Sc % global parameters, most of them are dimensionless

global  dt FT FS_interval % control parameters

% Open files to save computed data, it is saved every FS_interval time
% steps


fid1 = fopen('./data/S','w');      % file for S: substrate

fid2 = fopen('./data/X','w');      % file for X: cell concentration 

fid3 = fopen('./data/QS','w');  % file for QS: Quorum sensing molecues

fid4 = fopen('./data/WL','w');  % file for WL : lsaI

fid5 = fopen('./data/WR','w');  % file for WR : rsaL

y0 = Initial_QS();

sol = ode45(@QS_rhs, [0 FT], y0);  % solve it using ode45

dst = dt*FS_interval;

t_int = 0:dst:FT;

s_int = deval(sol,t_int);


fprintf(fid1, '%12.10e ', s_int(1,:));
fprintf(fid1, '\n');

fprintf(fid2, '%12.10e ', s_int(2,:));
fprintf(fid2, '\n');

fprintf(fid3, '%12.10e ', s_int(3,:));
fprintf(fid3, '\n');

fprintf(fid4, '%12.10e ', s_int(4,:));
fprintf(fid4, '\n');

fprintf(fid5, '%12.10e ', s_int(5,:));
fprintf(fid5, '\n');

fclose('all');
    
    