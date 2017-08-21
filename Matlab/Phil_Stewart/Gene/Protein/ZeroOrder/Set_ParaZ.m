% Script file to define the global variables, nondimensionalize the
% equations (by scaling the coefficients in the differentiatoin equations)
% Here the values are hand-typing, it should be obtained by reading
% external data file in the future
% Run the simulation and save the results

fprintf(1, '\n Start running simulation \n');

global  Mu_max Rho Dw Ds Sigma Kg LSc Sc % global parameters, most of them are dimensionless

d_s = 1.53e-9; % dimensional diffusion coefficient for nutrient (Oxygen)

d_w = 1e-16; % dimensional diffusion coefficient for GFP (can be set to zero if non diffusion is allowed)

k = 15;  % unit = mg/l/s coefficient in the zero-th order reaction kinetic for substrate consumption

mu_max = 0.1/3600; % maximum growth rate for bioflilm growth

rho = 6 ;  % unit = g/l total protein density

s0 = 6;  % unit = mg/l reference value substrate concentration

fi0 = 1 ;    % reference value for inducer concentration

wi0 = 6;    % unit = g/l reference value for protein density

kg = 0.5/3600; %0.053/3600; % unit = 1/h, decay rate of protein due to 

h0 = 1e-4;  % reference length scale : 100 micrometer

t0 = 3600;  % reference time scale : 1 hour

sigma = 0; % 0.01; % dimensional detachment coefficient

Dw = d_w*t0/h0/h0;

Ds = d_s*s0/k/h0/h0;

Mu_max = t0*mu_max;

Rho = rho*fi0/wi0;

Kg = kg*t0;

LSc = 0.8; % zeta < LSc : no substrate consumption, zeta >= LSc : with substrate consumption,  zero-th order kinetics

Sc = 1e-5; % S < Sc : no protein production, s >= Sc: with protein production , zero-th order kinetics

Sigma = sigma*h0*t0;

global Nx dx dt FT FS_interval flag_w flag_kg  % control parameters


Nx = 1000;  % number of subtervals in the mesh

dx = 1.0/Nx; % spatial grid size

dt = 0.005;   % temporal grid size

FT = 24; % final time of the simulation, it's 24 hours since t0 = 1 hour

FS_interval = FT/dt/24; % Save 100 intermediate steps

flag_w = 2; % 1 : with diffusion;  2: without diffusion

flag_kg = 0; % 1: with protein decay, 0: without protein decay