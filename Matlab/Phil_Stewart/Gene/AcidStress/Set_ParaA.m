% Script file to define the global variables, nondimensionalize the
% equations (by scaling the coefficients in the differentiatoin equations)
% Here the values are hand-typing, it should be obtained by reading
% external data file in the future
% Run the simulation and save the results

fprintf(1, '\n Start running simulation \n');

global  mu_max Rho Sigma Kg % global parameters, most of them are dimensionless

global s0 p0  h0 t0 PD G_s G_p

d_s = 2.24e-10; % dimensional diffusion coefficient for nutrient (glucose)

d_p = 3.51e-10; % dimensional diffusion coefficient for product (lactate)

mu_max = 0.8/3600; % maximum growth rate for bioflilm growth

s0 = 800;  % unit = mg/l reference value of glucose concentration, aslo the bulk glucose concentration

p0 = 500;  % unit = mg/l, reference value of lactate concentration

kg = 0.1/3600; %0.053/3600; % unit = 1/h, decay rate of protein due to 

h0 = 1e-4;  % reference length scale : 100 micrometer

t0 = 3600;  % reference time scale : 1 hour

PD = 1e-4;   % penetration depth of the nutrient

sigma = 0.08; % 0.01; % dimensional detachment coefficient

Y_ps = 0.9; % Yield coefficient of lactate on glucose: 0.9 mg lactate per mg glucose

Y_Xs = 0.5; %  Yield coefficient of biomass on glucose: 0.5 mg biomass per mg glucose

Xd = 50000;  % density of biomass, mg of X /l

k = 90;  % unit = mg/l/s coefficient in the zero-th order reaction kinetic for substrate consumption

G_s = k*(h0^2*Xd)/(d_s*s0*Y_Xs);  % reaction factor of glucose, G_s*mu is the reaction term, use k to adjust reaction rate

G_p = -k*(h0^2*Xd*Y_ps)/(d_p*p0*Y_Xs); % reaction factor of lactate, G_p*mu is the reaction term

rho = 6 ;  % unit = mg/l  gene density
wi0 = 6;   % unit = mg/l reference value for gene density
fi0 = 1;   % reference value for gene growth factor, dimensionless

Rho = rho*fi0/wi0;

Kg = kg*t0;

Sigma = sigma*h0*t0; % dimensionless detachment coefficient

global Nx dx dt FT FS_interval flag_kg  % control parameters


Nx = 200;  % number of subtervals in the mesh

dx = 1.0/Nx; % spatial grid size

dt = 0.005;   % temporal grid size

FT = 48; %24; % final time of the simulation, it's 24 hours since t0 = 1 hour

FS_interval = FT/dt/48; % Save 100 intermediate steps

flag_kg = 1; % 1: with protein decay, 0: without protein decay