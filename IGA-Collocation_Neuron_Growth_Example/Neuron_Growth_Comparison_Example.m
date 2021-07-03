% IGA-collocation Implementation for 2D neuron growth
% Kuanren Qian
% 06/02/2021

%% CleanUp
close all;
clear;
clc;

%% Including Path
addpath('./IGA_collocation_algorithm');

disp('********************************************************************');
disp('2D Phase-field Neuron Growth solver using IGA-Collocation');
disp('********************************************************************');

%% Variable Initialization
% time stepping variables
dtime = 1e-3; %actual dt is calculated in NR
end_iter = 10000;

% tolerance for NR method
tol = 1e-4;

% B-spline curve order (U,V direction)
p = 3;
q = 3;
Nx = 40;
Ny = 40;
dx = 1;
dy = 1;
knotvectorU = [0,0,0,linspace(0,Nx,81),Nx,Nx,Nx].';
knotvectorV = [0,0,0,linspace(0,Ny,81),Ny,Ny,Ny].';

% setting lenu lenv this way for easier access to ghost nodes later on
lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(p-1);

% neuron growth variables
abar = 0.45;
% ab = (abar/(1+delta));
aniso = 6;
k_conc = 0.5;
kappa= 4;
alph = 0.9; % changing name to alph cause alpha is a function
pix=4.0*atan(1.0);
gamma = 15.0;
tau = 0.005;
% tau = 0.0003;
delta_w = 4*dx;
gamma_w = 1;
lambda = 1.0e-16;
b = 2*atanh(1-2*lambda); %88.9
W = (6*gamma_w*b)/delta_w;
M = 10;
% M_phi = (sqrt(2*W)/(6*abar))*M; %50
M_phi = 1/tau;
M_theta = 0.5*M_phi;
D_cell = 10;
D_liquid = 1;
s_coeff = 0.002;
H = 1e-6;

delta = 0.2;
epsilonb = 0.05;

% Seed size
seed_radius = 10;

% Expanding domain parameters
BC_tol = 10;
expd_coef = 1.2;

disp('Base variable - initialization done!');

%% Iterating variable initialization
% initializing phi and concentration based on neuron seed position
[phi, conc, conct] = initialize_neurite_growth(seed_radius, lenu, lenv, dx);

% reshpae phi and concentration for calculation
phi = reshape(phi,lenu*lenv,1);
conc = reshape(conc,lenu*lenv,1);
conct = reshape(conct,lenu*lenv,1);

%% Constructing coef matrix
order_deriv = 2;    % highest order of derivatives to calculate
sprs = 1;   % sparse or not (for kqCollocationDers)
[NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points,dersU] ...
    = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
lap = N2uNv + NuN2v;

phi = NuNv\phi;
conc = NuNv\conc;
conct = NuNv\conct;

% initializing theta and temperature
theta=rand(lenu,lenv);
theta=reshape(theta,lenu*lenv,1);
theta = NuNv\theta;
tempr = zeros([lenu*lenv,1]);

phi_initial = reshape(phi,lenu,lenv);
conc_initial = reshape(conc,lenu,lenv);
theta_initial = reshape(theta,lenu,lenv);
tempr_initial = reshape(tempr,lenu,lenv);
for i = 2:lenu-1
    for j = 2:lenv-1
        phi_initial(i,j) = 0;
        conc_initial(i,j) = 0;
        theta_initial(i,j) = 0;
        tempr_initial(i,j) = 0;
    end
end
phi_initial = reshape(phi_initial,lenu*lenv,1);
conc_initial = reshape(conc_initial,lenu*lenv,1);
theta_initial  = reshape(theta_initial,lenu*lenv,1);
tempr_initial  = reshape(tempr_initial,lenu*lenv,1);

% plotting initial phi
set(gcf,'position',[700,100,700,900]);
% set(gcf,'position',[200,100,700,300]);
colormap parula;

% ID for boundary location (suppress 4 edges)
% id = 1 means there is bc
bcid = zeros([lenu,lenv]);
for i = 1:lenu
    bcid(1,i) = 1;
    bcid(lenu,i) = 1;
    bcid(i,1) = 1;
    bcid(i,lenv) = 1;
end
bcid = reshape(bcid,lenu*lenv,1);

disp('Phi,conc,theta,tempr,bcid - initialization done!');
disp('********************************************************************');

%% Transient iteration computation
disp('Starting Neuron Growth Model transient iterations...');

phi = sparse(phi);
conc = sparse(conc);
theta = sparse(theta);
tempr = sparse(tempr);
conct = sparse(conct);
phi_ones = sparse(zeros([lenu*lenv,1]));
bcid = sparse(bcid);

phi_initial_temp = reshape(phi,lenu,lenv);
[phidyo,phidxo] = gradient_mat(phi_initial_temp,Nx,Ny,dx,dy);
sq_grad_phi = sparse(reshape(phidyo.^2+phidxo.^2,lenu*lenv,1));
sum_grad_phi_ori = sum(sq_grad_phi);

% transient iterations
for iter=1:1:end_iter
    tic
    fprintf('Progress: %.2d/%.2d\n',iter,end_iter);

    % calculating a and a*a' (aap) in the equation using theta and phi
    [a, ~, aap,pdy,pdx] = kqGetEpsilonAndAap(epsilonb,delta,phi,theta,NuNv,NuN1v,N1uNv);
    a = reshape(a,lenu*lenv,1);
    aap = reshape(aap,lenu*lenv,1);

    teq = 1;
    E = (alph./pix).*atan(gamma.*(teq-(NuNv*tempr)));       

    %% Phi (Implicit Nonlinear NR method)
    % NR method initial guess (guess current phi)
    phiK = phi;
    % residual for NR method (make sure it is larger than tolerance at the beginning)
    R = 2*tol;
    
    % magnitude of theta gradient (offset by 1e-12 to avoid division by zero)
    mag_grad_theta = sparse(sqrt((N1uNv*theta).^2+(NuN1v*theta).^2)+1e-12);

    % NR method calculation
    ind_check = 0;
    NNa = NuNv*a;
    N1Na = N1uNv*a;
    NN1a = NuN1v*a;
    NNaap = NuNv*aap;
    N1Naap = N1uNv*aap;
    NN1aap = NuN1v*aap;

    dt_t = 0;
    while max(abs(R)) >= tol        
        NNpk = NuNv*phiK;
        N1Npk = N1uNv*phiK;
        NN1pk = NuN1v*phiK;   
        N1N1pk = N1uN1v*phiK;   
        LAPpk = lap*phiK;

        terma2 = 2*NNa.*N1Na.*N1Npk+NNa.^2.*LAPpk ...
            +2*NNa.*NN1a.*NN1pk;
        termadx = N1Naap.*NN1pk+NNaap.*N1N1pk;
        termady = NN1aap.*N1Npk+NNaap.*N1N1pk;     
        termNL = -NNpk.^3+(1.5-E).*NNpk.^2+(E-0.5).*NNpk;
%         termP1theta = (NNpk-NNpk.^2).*6.*H.*mag_grad_theta;
% 
        terma2_deriv =  2*NNa.*N1Na.*N1uNv+NNa.^2.*lap ...
            + 2*NNa.*NN1a.*NuN1v;
        %termadx_deriv = N1Naap.*NuN1v+NNaap.*N1uN1v;
        %termady_deriv = NN1aap.*N1uNv+NNaap.*N1uN1v;  
        termadxdy_deriv = N1Naap.*(-NuN1v+N1uNv);
        termNL_deriv = -3*NNpk.^2+2*(1.5-E).*NNpk+(E-0.5);
        termNL_deriv = termNL_deriv.*NuNv;
%         termP1theta_deriv = (1-2.*NNpk).*6.*H.*mag_grad_theta;
%         termP1theta_deriv = termP1theta_deriv.*NuNv;
%         
%         R = 1/tau*(terma2-termadx+termady+termNL-termP1theta);
        R = 1/tau*(terma2-termadx+termady+termNL);
        R = R*dtime-NNpk+(NuNv*phi);
%         dR = 1/tau*(terma2_deriv-termadx_deriv+termady_deriv+termNL_deriv-termP1theta_deriv);
        dR = 1/tau*(terma2_deriv+termadxdy_deriv+termNL_deriv);
        dR = dR*dtime-NuNv;

        % check residual and update guess
        R = R - dR*phi_initial;
        [dR, R] = StiffMatSetupBCID(dR, R,bcid,phi_initial);
        dp = dR\(-R);
        phiK = phiK + dp;
        
        max_phi_R = full(max(abs(R)));
        fprintf('Phi NR Iter: %.2d -> max residual: %.2d\n',ind_check, max_phi_R);
        if (ind_check >= 100 || max(abs(R))>1e20)
            error('Phi NR method NOT converging!-Max residual: %.2d\n',max_phi_R);
        end
        ind_check = ind_check + 1;
        dt_t = dt_t+dtime;
    end
 
    dtime_theta = dt_t/10;
    dtime_conc = dt_t;
    dtime_tempr = dt_t;
    dtime_conct = dt_t;

    %% Temperature (Explicit method)
    temprLHS = NuNv;
    temprRHS = (NuNv*tempr + lap*tempr.*dtime_tempr + kappa*(NuNv*phiK-NuNv*phi));
    temprRHS = temprRHS - temprLHS*tempr_initial;
    [temprLHS, temprRHS] = StiffMatSetupBCID(temprLHS, temprRHS,bcid,tempr_initial);
    tempr_new = temprLHS\temprRHS;
    
    %% Theta (Implicit method)
    lap_theta = lap*theta;
    P2 = 10*NNpk.^3-15*NNpk.^4+6*NNpk.^5;
    P2 = P2.*s_coeff.*mag_grad_theta;

    thetaLHS = (NuNv-dtime_theta.*M_theta.*P2.*lap);
    thetaRHS = (NuNv*theta);
    thetaRHS = thetaRHS - thetaLHS*theta_initial;    
    [thetaLHS, thetaRHS] = StiffMatSetupBCID(thetaLHS, thetaRHS,bcid,theta_initial);
    theta_new = thetaLHS\thetaRHS;

    %% iteration update
    % update variables in this iteration
    phi = phiK;
    theta = theta_new;
    tempr = tempr_new;
    
    %% Plotting figures
    if(mod(iter,10) == 0 || iter == 1)
        phi_plot = reshape(NuNv*phiK,lenu,lenv);
        subplot(4,2,1);
        imagesc(phi_plot);
        title(sprintf('Phi at iteration = %.2d',iter));
        axis square;
        colorbar;

        subplot(4,2,2);
        E_plot = reshape(NuNv*E,lenu,lenv);
        imagesc(E_plot);
        title(sprintf('E at iteration = %.2d',iter));
        axis square;
        colorbar;

        tempr_plot = reshape(NuNv*tempr_new,lenu,lenv);
        subplot(4,2,3);
        imagesc(tempr_plot(2:end-1,2:end-1));
        title(sprintf('Tempr at iteration = %.2d',iter));
        axis square;
        colorbar;

        theta_plot = reshape(NuNv*theta_new,lenu,lenv);
        subplot(4,2,4);
        imagesc(theta_plot);
        title(sprintf('theta_plot at iteration = %.2d',iter));
        axis square;
        colorbar;

        a_plot = reshape(NuNv*a,lenu,lenv);
        subplot(4,2,5);
        imagesc(a_plot);
        title(sprintf('a at iteration = %.2d',iter));
        axis square;
        colorbar;
      
        pdy_plot = reshape(pdy,lenu,lenv);
        pdx_plot = reshape(pdx,lenu,lenv);

        subplot(4,2,6);
        imagesc(pdy_plot);
        title(sprintf('pdy at iteration = %.2d',iter));
        axis square;
        colorbar;
         
        subplot(4,2,7);
        imagesc(pdx_plot);
        title(sprintf('pdx at iteration = %.2d',iter));
        axis square;
        colorbar;
        
        atan2_plot = reshape(atan2(pdy,pdx),lenu,lenv);
        subplot(4,2,8);
        imagesc(atan2_plot);
        title(sprintf('atan2 at iteration = %.2d',iter));
        axis square;
        colorbar;

        % plot current iteration
        drawnow;

        if(mod(iter,500) == 0)
            saveas(gcf,sprintf('NeuronGrowth_ex1_%.2d.png',iter));
        end
        
        %% This portion of code expands domain, not needed for comparison case
%         if(iter~=1)
%             if( max(max(phi_plot(1:BC_tol,:))) > 0.5 || ...
%                 max(max(phi_plot(:,1:BC_tol))) > 0.5 || ...
%                 max(max(phi_plot(end-BC_tol:end,:))) > 0.5 || ...
%                 max(max(phi_plot(:,end-BC_tol:end))) > 0.5)
%            
%                 disp('********************************************************************');
%                 disp('Expanding Domain...');
%                 
%                 Nx = floor(Nx*expd_coef);
%                 Ny = floor(Ny*expd_coef);
%                 
%                 knotvectorU = [0,0,0,linspace(0,Nx,Nx+1),Nx,Nx,Nx].';
%                 knotvectorV = [0,0,0,linspace(0,Ny,Ny+1),Ny,Ny,Ny].';
% 
%                 lenu = length(knotvectorU)-2*(p-1);
%                 lenv = length(knotvectorV)-2*(p-1);
% 
%                 oldNuNv = NuNv;
%                 [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points] ...
%                     = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs);
%                 lap = N2uNv + NuN2v;
% 
%                 sz = length(lap);
%                 [phi,theta,~,tempr,phi_initial,theta_initial,~,tempr_initial,bcid] ...
%                     = kqExpandDomain_Comparison_Example(sz,phiK,theta_new,theta_new,tempr_new,oldNuNv,NuNv);
%                 
%                 toc
%                 disp('********************************************************************');
%             end
%         end
    end
    toc
end

disp('All simulations complete!\n');