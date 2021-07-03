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

addpath('./utils');

%== get initial wall time:
time0 = clock();
format long;
%-- Simulation cell parameters:
Nx =80;
Ny =80;
NxNy = Nx*Ny;
% dx = 0.03;
% dy = 0.03;
dx = 0.5;
dy = 0.5;
ii = 1000;
%--- Time integration parameters:
nstep = 200000;
nprint = 100;
dtime = 1e-4;
flag = 1;
%--- Material specific parameters:
tau = 0.005;
% tau = 0.001;
epsilonb = 0.05;
mu = 1.0;
kappa = 4.0;
delta = 0.2;
aniso = 6.0;
alpha = 0.9;
gamma = 15.0;
teq = 1.0;
%theta0 = 3*(rand(Nx,Ny)-0.5);
xtheta = -pi/2;
% for i=1:Nx
%     for j=1:Ny
%         incl = atan2((i-Ny/2),(j-Nx/2));
%         %theta0(i,j) = 0.5 + 0.5*incl/(pi);
%         if(incl>xtheta-pi/12 && incl <xtheta+pi/12)
%         elseif(incl>0-pi/12 && incl <0+pi/12)    
%         else
%             incl = 3*(rand(1,1)-0.5);
%         end
%         %    theta1(i,j) = 0.5 + 0.5*incl/(pi);
%         %else
%         theta0(i,j) = incl;
%         %end
%     end
% end
theta0 = 2*pi*rand(Nx);

% seed = 1100;
seed = 10;
seed = seed^2;
%
pix=4.0*atan(1.0);
%--- Initialize and introduce
% initial nuclei:
[phi,tempr] = nucleus(Nx,Ny,seed);

%---
%--- Evolution
%---
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
    %-- calculate angle:
%     theta = atan2(phidy,phidx)+pi/2;
    theta = atan2(phidy,phidx);
    %--- epsilon and its derivative:
    epsilon = epsilonb*(1.0+delta*cos(aniso*(theta-theta0)));
    epsilon_deriv = -epsilonb*aniso*delta*sin(aniso.*(theta-theta0));
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
        phiold.*(1.0-phiold).*(phiold - 0.5 + m));
    
    %-- evolve temperature:
    tempr = tempr + dtime*lap_tempr + kappa*(phi-phiold);
    
%     if(mod(istep,500)==0)
%     for i =1:Nx
%         for j=1:Ny
%             incl = atan2((i-Ny/2),(j-Nx/2));
%             if(incl>xtheta-pi/12 && incl <xtheta+pi/12)
%                 theta0(i,j) = theta0(i,j) + (-1)^flag*pi/12;
%             end
%             if(incl>0-pi/12 && incl <0+pi/12)
%                 theta0(i,j) = theta0(i,j) + (-1)^flag*pi/12;
%             end
%         end
%     end
%     end
%     if(mod(istep,ii)==0)
%         flag = flag+1;
%         ii = ii+2000;
%     end
    %---- print results

    if(mod(istep,nprint) == 0 )
        saveas(gcf,sprintf('NeuronGrowth_ex1_%.2d.png',istep));
        fprintf('done step: %5d\n',istep);
        subplot(2,2,1)
        imagesc(phi)
        title(sprintf('phi at iteration = %.2d',istep))
        colorbar
        subplot(2,2,2)
        imagesc(tempr)
        title("tempr")
        colorbar
        subplot(2,2,3)
        imagesc(theta0)
        title("theta")
        colorbar
        drawnow 
    end %if
end %istep
%--- calculate compute time:
compute_time = etime(clock(), time0);
fprintf('Compute Time: %10d\n', compute_time);



