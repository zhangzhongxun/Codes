% Script file to define the global variables, nondimensionalize the
% equations (by scaling the coefficients in the differentiatoin equations)
% Here the values are hand-typing, it should be obtained by reading
% external data file in the future
% Run the simulation and save the results

fprintf(1, '\n Start running simulation \n');

global  Mu_max  Dq Ds Sigma Kg_WL Kg_WR Kg_Q LSc Sc t0 K_s K_ex % global parameters, most of them are dimensionless

d_s = 1.53e-9; % dimensional diffusion coefficient for nutrient (Oxygen)

d_q = 1.2e-10; % dimensional diffusion coefficient for QS molecules)

k_s = 15;  % unit = mg/l/s coefficient in the zero-th order reaction kinetic for substrate consumption

mu_max = 0.2/3600; % maximum growth rate for bioflilm growth

rho = 6 ;  % unit = mg/l total protein density

s0 = 6;  % unit = mg/l reference value substrate concentration

fi0 = 1 ;    % reference value for inducer concentration

%wi0 = 6;    % unit = mg/l reference value for protein density

kg_WL = 0.053/3600; % unit = 1/h, turnover rate of lsaI

kg_WR = 0.053/3600; % unit = 1/h, turnover rate of rsaL

kg_Q = 0.053/3600; % unit = 1/h, turnover rate of QS

h0 = 1e-4;  % reference length scale : 100 micrometer

t0 = 3600.0;  % reference time scale : 1 hour

sigma = 0; % 0.01; % dimensional detachment coefficient

K_s = k_s*t0;

k_ex = 0.01;  % exchange rate of oxygen between batch and air

K_ex = k_ex*s0;

Dq = d_q*t0/h0/h0;

Ds = d_s*s0/k_s/h0/h0;

Mu_max = t0*mu_max;

Kg_WL = kg_WL*t0;

Kg_WR = kg_WR*t0;

Kg_Q = kg_Q*t0;

LSc = 0.8; % zeta < LSc : no substrate consumption, zeta >= LSc : with substrate consumption,  zero-th order kinetics

Sc = 0e-7; % S < Sc : no protein production, s >= Sc: with protein production , zero-th order kinetics

Sigma = sigma*h0*t0;

global QS_0 WL_0 WR_0 Alpha_1 Alpha_2 Alpha_3 Alpha_4 max_rate_WL Q_bar Q_BC g_l g_r g_x

X0 = 0.5;
WL_0 = 6.0;
WR_0 = 6.0;
QS_0 = 6.0;
rho = 60.0; % protein density

alpha_1 = 6e4;
alpha_2 = 1e2; 
alpha_3 = 3;
alpha_4 = 1;

Alpha_1 = alpha_1*X0*WL_0/QS_0;
Alpha_2 = alpha_2*WR_0;
Alpha_3 = alpha_3*QS_0/WR_0;
Alpha_4 = alpha_4*s0;

Y_xs = 0.5;
g_x = X0/Y_xs; 

g_l = rho/WL_0;
g_r = rho/WR_0;

max_rate_WL = 1;  %10/3600;
Q_bar = 0.1;
Q_BC = 0; % value of Dirichlet BC at x = 1

global Nx dx dt FT FS_interval flag_kg_WL  flag_kg_WR flag_kg_Q % control parameters


Nx = 200;  % number of subtervals in the mesh

dx = 1.0/Nx; % spatial grid size

dt = 0.02;   % temporal grid size

FT = 72; % 72; % final time of the simulation, it's 24 hours since t0 = 1 hour

FS_interval = FT/dt/72; % Save 100 intermediate steps

flag_kg_Q = 0; % 1: with QS turnover, 0: without QS turnover

flag_kg_WL = 0; % 1: with lsaI turnover, 0: without lsaI turnover

flag_kg_WR = 0; % 1: with with rsaL turnover, 0: without rsaL turnover