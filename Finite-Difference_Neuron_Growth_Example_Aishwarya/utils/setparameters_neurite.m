function parameters = setparameters_neurite()

maxlevel = 3;

rho = zeros(maxlevel,1);
nstep = zeros(maxlevel,1);

rho(1,1) = 1;
rho(2,1) = 1;
rho(3,1) = 3;

nstep(1,1) = 600;
nstep(2,1) = 900;
nstep(3,1) = 500;
nprint = 50;

nelemx = 80;
nelemy = 80;

pU = 3;
pV = 3;

Nx = 400;
Ny = 400;

dx = 1;
dy = 1;

%% Time integration parameters:
dtime = 5e-3;

%% Material specific parameters:
kappa= 1.8;
delta=0.5;
aniso= 2.0;
alpha = 0.9;
gamma = 15.0;
seed = (40*dx)^2;
k_conc = 0.5;
D_cell = 10;
D_liquid = 1;

%% Calculate constant values

abar = 0.45;
epsilonb = (abar/(1+delta));
tau = epsilonb;
M_phi =50;
M_theta = 0.5*M_phi ;

%% Tubulin parameters
alpha_t = 0.001;
beta_t = 0.001;
Diff = 2.0;
source_coeff = 1.0;
s_coeff = 3.e-6;

flag_replicate = 0;

parameters = struct('pU',pU,'pV',pV,'maxlevel',maxlevel,'nelemx',nelemx,'nelemy',nelemy,'rho',rho,'Nx',Nx,'Ny',Ny,'dx',dx,'dy',dy,...
    'dtime',dtime,'kappa',kappa,'delta',delta,'aniso',aniso,'alpha',alpha,'gamma',gamma,'seed',seed,'k_conc',k_conc,'D_cell',D_cell,...
    'D_liquid', D_liquid,'abar',abar,'epsilonb',epsilonb,'tau',tau,'M_phi',M_phi,'M_theta', M_theta,'alpha_t',alpha_t,'beta_t',beta_t,...
    'Diff',Diff,'source_coeff',source_coeff,'flag_replicate',flag_replicate,'nstep',nstep,'s_coeff',s_coeff, 'nprint', nprint);
end