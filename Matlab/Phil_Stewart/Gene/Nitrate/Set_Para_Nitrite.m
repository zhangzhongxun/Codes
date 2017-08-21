% Script file to define the global variables, nondimensionalize the
% equations (by scaling the coefficients in the differentiatoin equations)
% Here the values are hand-typing, it should be obtained by reading
% external data file in the future
% Run the simulation and save the results

fprintf(1, '\n Start running simulation \n');

global  Mu_max_O2 Mu_max_NO2 Alpha D_O2 D_NO2 Kg  % global parameters, most of them are dimensionless

d_o2 = 1.53e-9; % dimensional diffusion coefficient for Oxygen in biofilm, unit = m^2/s
d_no2 = 1.32e-9; % dimensional diffusion coefficient for Oxygen in biofilm, unit = m^2/s

%d_w = 1e-16; % dimensional diffusion coefficient for GFP (can be set to zero if non diffusion is allowed)

k_o2 = 4/3.92;  % unit = 1/s, coefficient in the first order reaction kinetic for oxygen consumption
                % adjusted such that SO2(2) = 1/175

k_no2 = 4; % unit = 1/s, coefficient in the first order reaction kinetic for nitrite consumption

mu_max_o2 = 0.1/3600; % maximum growth rate for bioflilm growth when consume oxygen

mu_max_no2 = 0.1/3600; % maximum growth rate for bioflilm growth when consume nitrite

kg = 0.053/3600; % unit = 1/h, decay rate of RNA due to 

alpha = 6;  % protein density

%s0 = 1e-3;  % reference value substrate concentration

fi0 = 2 ;    % reference value for inducer concentration

wi0 = 6;    % unit = mg/l reference value for protein density

h0 = 1e-4;  % reference length scale : 100 micrometer

t0 = 3600;  % reference time scale : 1 hour

D_O2 = d_o2/k_o2/h0/h0;

D_NO2 = d_no2/k_no2/h0/h0;

Mu_max_O2 = t0*mu_max_o2;

Mu_max_NO2 = t0*mu_max_no2;

Alpha = alpha*fi0/wi0;

Kg = kg*t0;

global SO2_crit eta; % control parameters for RNA producing kinetics

SO2_crit = 1/175;

eta = log(2)/SO2_crit;

global Nx dx dt FT FS_interval flag_w flag_kg x1_O x2_O x1_N x2_N x3_N % control parameters

x1_O = 4;  % spatial coordinate where BCs for Oxygen concentration are specified 
x2_O = 0;

x1_N = 4;  % spatial coordinate where BCs for Nitrite concentration are specified 
x2_N = 2;
x3_N = 0;

Nx = 400;  % number of subtervals in the mesh

dx = x1_O/Nx; % spatial grid size

dt = 0.01;   % temporal grid size

FT = 48; % final time of the simulation, it's 100 days since t0 = 1 hour

FS_interval = FT/dt/48; % Save 24 intermediate steps

flag_w = 2; % 1 : with diffusion;  2: without diffusion

flag_kg = 0; % 0 : no RNA turnover, 1: with RNA turnover
