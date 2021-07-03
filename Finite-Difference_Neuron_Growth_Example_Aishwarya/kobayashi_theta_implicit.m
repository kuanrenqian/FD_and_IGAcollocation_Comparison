%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE-FIELD FINITE-DIFFRENCE %
%
% CODE FOR
%
% DENDRITIC SOLIDIFICATION
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%== get initial wall time:

clc
clear all
close all

time0 = clock();
format long;
%-- Simulation cell parameters:
Nx = 100;
Ny = 100;
NxNy = Nx*Ny;
dx = 0.03;
dy = 0.03;
%--- Time integration parameters:
nstep = 6000;
nprint = 2;
dtime = 1e-4;

%--- Material specific parameters:
tau = 0.0003;
epsilonb = 0.01;
mu = 1.0;
kappa = 1.8;
delta = 0.02;
aniso = 6.0;
alpha = 0.9;
gamma = 10.0;
teq = 1.0;
theta00 = 0.2;
seed = 50.0;
gamma_chi = 800;
s_coeff1 = 0.002;
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = 1.0;
H = 1e-6;
M_phi = 1/tau;
M_theta = 0.5*M_phi;

pix=4.0*atan(1.0);
%--- Initialize and introduce
% initial nuclei:
[phi,tempr,conc_t] = nucleus(Nx,Ny,seed);
theta = rand(Nx,Ny);
theta0 = theta;
%--
%--- Evolution
%---
[phidyo,phidxo] = gradient_mat_imf(phi,dx,dy);
sq_grad_phi = phidyo.^2+phidxo.^2;
sum_grad_phi_ori = sum(sq_grad_phi(:));

for istep =1:nstep
    phiold =phi;
    %---
    % calculate the laplacians
    %and epsilon:
    %---
    lap_phi = laplacian_imf(phi, dx, dy);
    %--
    lap_tempr = laplacian_imf(tempr, dx, dy);
    %--gradients of phi:
    [phidy,phidx]=gradient_mat_imf(phi,dx,dy);
    [thetady,thetadx]=gradient_mat_imf(theta,dx,dy);
    theta_absolute = sqrt(thetadx.^2+thetady.^2);
    %-- calculate angle:
    atheta =atan2(phidy,phidx);
    xtheta = 2.*pi.*theta;
    %--- epsilon and its derivative:
    epsilon = epsilonb*(1.0+delta*cos(aniso*(atheta-theta00)));
    epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(atheta-theta00));
    %--- first term:
    dummyx =epsilon.*epsilon_deriv.*phidx;
    [term1,~] =gradient_mat_imf(dummyx,dx,dy);
    %--- second term:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [~,term2] =gradient_mat_imf(dummyy,dx,dy);
    %--- factor m:
    m =(alpha/pix)*atan(gamma*(teq-tempr));
    
    %-- Time integration:
    phi = phi +(dtime/tau) *(term1 +term2 + epsilon.^2 .* lap_phi + ...
        phiold.*(1.0-phiold).*(phiold - 0.5 + m))-(6.*phi.*(1-phi).*H.*theta_absolute);
    
    A = zeros(Nx*Ny, Nx*Ny);
    A(1,1) = 1;
    A(Nx*Ny, Nx*Ny) = 1;
    
    B = zeros(Nx*Ny, 1);
    B(1,1) = theta0(1,1);
    B(Nx*Ny,1) = theta0(1,1);
    
    [thetadx,thetady] = gradient_mat_imf(theta,dx,dy);
    [conctdy,conctdx]=gradient_mat_imf(conc_t,dx,dy);
    theta_absolute = sqrt(thetadx.^2 + thetady.^2);
    chi_gamma1 = zeros(size(theta));
    chi_gamma2 = zeros(size(theta));
    
%     %evaluate d_t*phi*gradient(conc_tubulin)
%     term_ct1x = phi.*conctdx;
%     term_ct1y = phi.*conctdy;
%     
%     %calculate divergence of d_t*phi*gradient(conc_tubulin)
%     [dummy,termctx] = gradient_mat_imf(term_ct1x,dx,dy);
%     [termcty,dummy] = gradient_mat_imf(term_ct1y,dx,dy);
%     
%     %evaluate gradient(phi*conc_tubulin)
%     phi_conc = phi.*conc_t;
%     [termpcy,termpcx] = gradient_mat_imf(phi_conc,dx,dy);
    
    
    for i =1:size(theta,1)
        for j =1:size(theta,2)
            %(u(i+1)-u(i))/dx
            %(chi_gamma(gamma, ((u(i+1)-u(i))/dx)))
            
            if(i==1)
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j) = 1;
                B((i-1)*(size(theta,1))+j,1) = theta0(1,1);
                
            elseif(j==1)
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j) = 1;
                B((i-1)*(size(theta,1))+j,1) = theta0(1,1);
                
            elseif(i==size(theta,1))
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j) = 1;
                B((i-1)*(size(theta,1))+j,1) = theta0(1,1);
                
            elseif(j==size(theta,2))
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j) = 1;
                B((i-1)*(size(theta,1))+j,1) = theta0(1,1);
            else
                p_ij = phi(i,j)^2*(3-2*phi(i,j));
                p_i1j = phi(i+1,j)^2*(3-2*phi(i+1,j));
                p_i2j = phi(i-1,j)^2*(3-2*phi(i-1,j));
                p_ij1 = phi(i,j+1)^2*(3-2*phi(i,j+1));
                p_ij2 = phi(i,j-1)^2*(3-2*phi(i,j-1));
                p_ihalfj = -M_theta*H*0.5*(p_ij+p_i1j);
                p_ihalf2j = -M_theta*H*0.5*(p_ij+p_i2j);
                p_ijhalf = -M_theta*H*0.5*(p_ij+p_ij1);
                p_ijhalf2 = -M_theta*H*0.5*(p_ij+p_ij2);
                chi_gamma1(i,j) = chi_gamma(gamma,theta_absolute(i,j));
                chi_gamma2(i,j) = chi_gamma(gamma,theta_absolute(i,j));
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j) = (1/dtime) + (p_ihalfj/dx^2)*chi_gamma1(i,j) + (p_ihalf2j/dx^2)*chi_gamma2(i,j) + (p_ijhalf/dy^2)*chi_gamma1(i,j) + (p_ijhalf2/dy^2)*chi_gamma2(i,j);
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j-1) = -(p_ijhalf2/dx^2)*chi_gamma1(i,j);
                A((i-1)*(size(theta,1))+j,(i-1)*(size(theta,1))+j+1) = -(p_ijhalf/dx^2)*chi_gamma1(i,j);
                A((i-1)*(size(theta,1))+j,(i-2)*(size(theta,1))+j) = -(p_ihalf2j/dx^2)*chi_gamma1(i,j);
                A((i-1)*(size(theta,1))+j,(i)*(size(theta,1))+j) = -(p_ihalfj/dx^2)*chi_gamma1(i,j);
                B((i-1)*(size(theta,1))+j,1) = theta(i,j)/dtime;
            end
        end
    end
    
    theta1 = A\B;
    theta = reshape(theta1,[size(theta,1),size(theta,2)]);

    
    %-- evolve temperature:
    tempr = tempr + dtime*lap_tempr + kappa*(phi-phiold);
    
%     %evolve concentration of tubulin
%     termct_diffusion = Diff.*(termctx(i,j)+termcty(i,j));
%     termct_activetr = alpha_t*(termpcx(i,j)+termpcy(i,j));
%     termct_degradation = beta_t*phi(i,j)*conc_t(i,j);
%     termct_source = source_coeff*sq_grad_phi(i,j)/sum_grad_phi_ori;
%     %termct_final = lap_conct(i,j) + lap_phi(i,j); % + alpha_t*phiold(i,j)*(conctdy(i,j)+conctdx(i,j)) + beta_t*phiold(i,j)*conc_t(i,j);
%     conc_t(i,j) = conc_t(i,j)  + ((dtime*1e-3))*(termct_diffusion-termct_degradation+termct_source-termct_activetr);
    %---- print results
    if(mod(istep,nprint) == 0 )
        fprintf('done step: %5d\n',istep);
        subplot(2,2,1)
        imagesc(phi)
        title("\phi")
        colorbar
        subplot(2,2,2)
        imagesc(tempr)
        title("tempr")
        colorbar
        subplot(2,2,3)
        imagesc(theta)
        title("theta")
        colorbar
        drawnow
    end %if
    if(mod(istep,500) == 0)       
        str1 = sprintf('phi%d.mat',istep);
        save(str1, 'phi');
    end
end %istep
%--- calculate compute time:
compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);
