%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PHASE-FIELD FINITE-DIFFRENCE %
%
% CODE FOR
%
% DENDRITIC SOLIDIFICATION
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%== get initial wall time:
time0 = clock();
format long;
%-- Simulation cell parameters:
Nx = 500;
Ny = 500;
NxNy = Nx*Ny;
dx = 0.03;
dy = 0.03;
%--- Time integration parameters:
nstep = 60000;
nprint = 600;
dtime = 1e-4;
dtime1 = 1e-5;
cc=1;
%--- Material specific parameters:
tau = 0.001;
Mphi = 1/tau;
M_theta = 0.5*Mphi;
epsilonb = 0.01;
mu = 1.0;
kappa = 5.0;
delta = 0.5;
aniso = 2.0;
alpha = 0.9;
gamma = 10.0;
teq = 1.0;
gamma1 = 800;

seed = 200.0;
s_coeff1 = 0.002;
alpha_t = 0.001;
beta_t = 0.001;
Diff = 1e2;
source_coeff = 1.0;
H = 1e-6;
therm_diff = 2;

pix=4.0*atan(1.0);
%--- Initialize and introduce
% initial nuclei:
[phi,tempr,theta,conc_t] = nucleus_theta(Nx,Ny,seed);

theta0 = theta;

figure
subplot(1,3,1)
imagesc(phi)
title("\phi")
colorbar
subplot(1,3,2)
imagesc(tempr)
title("tempr")
colorbar
subplot(1,3,3)
imagesc(theta)
title("theta")
colorbar
drawnow

[phidyo,phidxo] = gradient_mat_imf(phi,dx,dy);
sq_grad_phi = phidyo.^2+phidxo.^2;
sum_grad_phi_ori = sum(sq_grad_phi(:));
%---
%--- Evolution
%---
figure
for istep =1:nstep
    phiold =phi;
    theta_old = theta;
    %---
    % calculate the laplacians
    %and epsilon:
    %---
    lap_phi = laplacian_imf(phi, dx, dy);
    %--
    lap_tempr = laplacian_imf(tempr, dx, dy);
    %--
    lap_theta = laplacian_imf(theta_old, dx, dy);
    
    %--gradients of phi:
    [phidy,phidx]=gradient_mat_imf(phi,dx,dy);
    [conctdy,conctdx]=gradient_mat_imf(conc_t,dx,dy);
    [thetady,thetadx]=gradient_mat_imf(theta,dx,dy);
    theta_absolute = sqrt(thetadx.^2 + thetady.^2);
    
    %-- calculate angle:
    atheta = atan2(phidy,phidx)+pi/2;
    
    %if(istep<500)
    %xtheta = 0;
    %   xtheta = 2.*pi.*theta+pi/2;
    %     if(istep>5000)
    %         xtheta = pi.*theta;
    %     else
    for i = 1:size(phi,1)
        for j = 1:size(phi,2)
            %if(j<=Nx/2)
            %incl = atan2((i-Ny/2),(j-Nx/2))+pi/2;
            %xtheta(i,j) = incl + pi/2*theta(i,j);
            xtheta(i,j) = 2*pi*theta(i,j);
            %else
            %    xtheta(i,j) = 0;
            %end
            
        end
    end
    %     end
    %xtheta = 0;
    %xtheta = 2.*pi.*theta;
    
    %--- epsilon and its derivative:
    epsilon = epsilonb*(1.0+delta*cos(aniso*(atheta-xtheta)));
    epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(atheta-xtheta));
    %--- first term:
    dummyx =epsilon.*epsilon_deriv.*phidx;
    [term1,~] =gradient_mat_imf(dummyx,dx,dy);
    %--- second term:
    dummyy =-epsilon.*epsilon_deriv.*phidy;
    [~,term2] =gradient_mat_imf(dummyy,dx,dy);
    %--- factor m:
    [thetadx,thetady] = gradient_mat_imf(theta,dx,dy);
    %p(phi) calculation
    theta_absolute = sqrt(thetadx.^2+thetady.^2)+1e-6;
    p_phi = (phiold.^2).*(3-2.*phiold);
    %[g_p_phix, g_p_phiy] =  gradient_mat_imf(p_phi, dx, dy);
    
    dummyx = p_phi.*s_coeff1.*(thetadx./theta_absolute);
    [termx,~] = gradient_mat_imf(dummyx,dx,dy);
    
    dummyy = p_phi.*s_coeff1.*(thetady./theta_absolute);
    [~,termy] = gradient_mat_imf(dummyy,dx,dy);
    
    q_phi = phi.^2*(1-phi).^2;
    
    %evaluate d_t*phi*gradient(conc_tubulin)
    term_ct1x = phi.*conctdx;
    term_ct1y = phi.*conctdy;
    
    %calculate divergence of d_t*phi*gradient(conc_tubulin)
    [dummy,termctx] = gradient_mat_imf(term_ct1x,dx,dy);
    [termcty,dummy] = gradient_mat_imf(term_ct1y,dx,dy);
    
    %evaluate gradient(phi*conc_tubulin)
    phi_conc = phi.*conc_t;
    [termpcy,termpcx] = gradient_mat_imf(phi_conc,dx,dy);
    
    for i = 1:size(phi,1)
        for j = 1:size(phi,2)
            m =(alpha/pix)*atan(gamma*(teq-tempr(i,j)));
            
            phi(i,j) = phi(i,j) +(dtime/tau) *(term1(i,j) +term2(i,j) + epsilon(i,j).^2 .* lap_phi(i,j) + phiold(i,j).*(1.0-phiold(i,j)).*(phiold(i,j) - 0.5 + m   + (6*s_coeff1).*theta_absolute(i,j)));
            %-- evolve temperature:
            if(i==1 || i==Nx || j==1 || j==Ny)
                theta(i,j) = theta0(i,j);
            else
                if(phi(i,j)>0.5)
                    theta(i,j) = theta(i,j);
                else
                    p_ij = phi(i,j)^2*(3-2*phi(i,j));
                    p_i1j = phi(i+1,j)^2*(3-2*phi(i+1,j));
                    p_i2j = phi(i-1,j)^2*(3-2*phi(i-1,j));
                    p_ij1 = phi(i,j+1)^2*(3-2*phi(i,j+1));
                    p_ij2 = phi(i,j-1)^2*(3-2*phi(i,j-1));
                    p_ihalfj = 0.5*(p_ij+p_i1j);
                    p_ihalf2j = 0.5*(p_ij+p_i2j);
                    p_ijhalf = 0.5*(p_ij+p_ij1);
                    p_ijhalf2 = 0.5*(p_ij+p_ij2);
                    chi_gamma1 = chi_gamma(gamma1,theta_absolute(i,j));
                    chi_gamma2 = chi_gamma(gamma1,theta_absolute(i,j));
                    theta(i,j) = theta(i,j) + (-dtime1*M_theta*H)*(p_ihalfj*chi_gamma1*(theta(i+1,j)-theta(i,j))*(1/dx^2) - p_ihalf2j*chi_gamma1*(theta(i,j)-theta(i-1,j))*(1/dx^2) + p_ijhalf*chi_gamma1*(theta(i,j+1)-theta(i,j))*(1/dy^2) - p_ijhalf2*chi_gamma1*(theta(i,j)-theta(i,j-1))*(1/dy^2));
                end
            end
            
            tempr(i,j) = tempr(i,j) + dtime*therm_diff*lap_tempr(i,j) + kappa*(phi(i,j)-phiold(i,j));
            
            %evolve concentration of tubulin
            termct_diffusion = Diff.*(termctx(i,j)+termcty(i,j));
            termct_activetr = alpha_t*(termpcx(i,j)+termpcy(i,j));
            termct_degradation = beta_t*phi(i,j)*conc_t(i,j);
            termct_source = source_coeff*sq_grad_phi(i,j)/sum_grad_phi_ori;
            %termct_final = lap_conct(i,j) + lap_phi(i,j); % + alpha_t*phiold(i,j)*(conctdy(i,j)+conctdx(i,j)) + beta_t*phiold(i,j)*conc_t(i,j);
            conc_t(i,j) = conc_t(i,j)  + ((dtime*1e-3))*(termct_diffusion-termct_degradation+termct_source-termct_activetr);
        end
    end
    
    %---- print results
%     if(istep==800)
%         for i = 1:Nx
%             for j = 1:Ny
%                 incl = atan2((i-Ny/2),(j-Nx/2));
%                 if(i>=Nx/2)
%                     incl = -incl;
%                 end
%                 theta(i,j) = 0.5+ 0.5*incl/(pi);
%             end
%         end
%     end
    %if(mod(istep,nprint) == 0 )
    %         cc = cc*-1;
    %
    %         if(i>=Nx/2)
    %             theta = theta + cc*0.03*pi/4;
    %         else
    %             theta = theta - cc*0.03*pi/4;
    %         end
    %     end
    
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
    subplot(2,2,4)
    imagesc(conc_t)
    title("conc_t")
    colorbar
    drawnow
    
    if(mod(istep,500) == 0 )
        str1 = sprintf('phi%d.mat',istep);
        str2 = sprintf('conct%d.mat',istep);
        save(str1, 'phi');
        save(str2, 'conc_t');
    end %if
end %istep
%--- calculate compute time:
compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);



