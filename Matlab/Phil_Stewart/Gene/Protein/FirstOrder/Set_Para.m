% Script file to define the global variables, nondimensionalize the
% equations (by scaling the coefficients in the differentiatoin equations)
% Here the values are hand-typing, it should be obtained by reading
% external data file in the future
% Run the simulation and save the results

fprintf(1, '\n Start running simulation \n');

global  Mu_max Alpha Dw Ds Sigma  % global parameters, most of them are dimensionless

d_s = 1.53e-9; % dimensional diffusion coefficient for nutrient

d_w = 1e-16; % dimensional diffusion coefficient for GFP (can be set to zero if non diffusion is allowed)

k = 5;  % unit = 1/s, coefficient in the first order reaction kinetic for substrate consumption

mu_max = 0.2/3600; % maximum growth rate for bioflilm growth

alpha = 6;  % protein density

%s0 = 1e-3;  % reference value substrate concentration

fi0 = 1 ;    % reference value for inducer concentration

wi0 = 6;    % unit = mg/l reference value for protein density

h0 = 1e-4;  % reference length scale : 100 micrometer

t0 = 3600;  % reference time scale : 1 hour

sigma = 0; % 0.01; % dimensional detachment coefficient

Dw = d_w*t0/h0/h0;

Ds = d_s/k/h0/h0;

Mu_max = t0*mu_max;

Alpha = alpha*fi0/wi0;

Sigma = sigma*h0*t0;

global Nx dx dt FT FS_interval flag_w  % control parameters


Nx = 100;  % number of subtervals in the mesh

dx = 1.0/Nx; % spatial grid size

dt = 0.01;   % temporal grid size

FT = 24; % final time of the simulation, it's 100 days since t0 = 1 hour

FS_interval = FT/dt/24; % Save 24 intermediate steps

flag_w = 2; % 1 : with diffusion;  2: without diffusion
