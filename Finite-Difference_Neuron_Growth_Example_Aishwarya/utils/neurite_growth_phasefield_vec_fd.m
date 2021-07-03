function [phi,conc,conc_t,theta] = neurite_growth_phasefield_vec_fd(param,phi,conc,conc_t,tempr,theta)

Nx = param.Nx;
Ny = param.Ny;
dx = param.dx;
dy = param.dy;

iteration=0;

%% This script runs the iteration loop
tic
for istep = 1:sum(param.nstep)
    
    %for plotting the iteration versus N
    iteration = iteration+1;
    
    
    [phidx, phidy] = gradient_mat_imf(phi, dx,dy);
    
    %store the previous iteration phi matrix
    phiold =phi;
    
    
    %---
    % calculate the laplacians and epsilon:
    %---
    
    [lap_phi]=laplacian_imf(phi, dx, dy);
    
    %--
    [lap_conc]=laplacian_imf(conc,dx,dy);
    
    %--
    [lap_tempr]=laplacian_imf(tempr, dx, dy);
    num_theta_inner_steps =1;
    for theta_step =1:num_theta_inner_steps
        %--gradients of phi, conc, cont_tubulin and theta:
        [thetadx,thetady] = gradient_mat_imf(theta,dx,dy);
        
        %evaluate absolute value of theta
        theta_absolute = sqrt(thetadx.^2+thetady.^2)+1e-6;
        
        %-- calculate angle:
        atheta = atan2(phidy,phidx);
        xtheta = 2*pi*theta;
        
        epsilon = param.epsilonb*(1.0+param.delta*cos(param.aniso*(atheta-xtheta)));
        epsilon_deriv = -param.epsilonb*param.aniso*param.delta*sin(param.aniso.*(atheta-xtheta));
        
        %calculate first term in phi equation:
        dummyx = epsilon.*epsilon_deriv.*phidx;
        [term1,~] = gradient_mat_imf(dummyx,dx,dy);
        
        %calculate second term in phi equation:
        dummyy =-epsilon.*epsilon_deriv.*phidy;
        [~,term2] =gradient_mat_imf(dummyy,dx,dy);
        
        %p(phi) calculation
        p_phi = (phi.^3).*(10-15.*phi+6.*(phi.^2));
        
        %calculate terms for theta evaluation
        s_coeff1 = param.s_coeff.*ones(Nx,Ny);
        dummyx = p_phi.*s_coeff1.*(thetadx./theta_absolute);
        [termx,~] = gradient_mat_imf(dummyx,dx,dy);
        
        dummyy = p_phi.*s_coeff1.*(thetady./theta_absolute);
        [~,termy] = gradient_mat_imf(dummyy,dx,dy);
        
        %evolve theta
        theta = theta + param.M_theta*(param.dtime/num_theta_inner_steps)*(termx+termy);
        
    end
    
    %delta_L
    term_change= ones(Nx,Ny);
    
    %delta T
    teq = -1.*conc+2;
    
    %E term in PHI evolution
    E = (param.alpha/pi)*term_change.*atan(param.gamma*(teq-tempr));
    
    % calculate q(phi)
    q_phi = phi.^2*(1-phi).^2;
    
    %evolve phi:
    phi = phi + param.M_phi*(param.dtime/param.tau)*(term1 +term2 + ...
        epsilon.^2.*lap_phi + phiold.*(1.0-phiold).*(phiold - 0.5 + E)-...
        (30*s_coeff1).*q_phi.*theta_absolute);
    
    %evolve temperature:
    tempr =tempr + param.M_phi*param.dtime*lap_tempr + param.kappa*(phi-phiold);
    
    %evolve concentration
    D = param.D_cell+(param.D_liquid-param.D_cell)*(1-phi)./(1-phi+param.k_conc*phi);
    termc = ((1-param.k_conc)*conc)./(1-phi+param.k_conc*phi);
    termc_final = D.*(lap_conc+termc.*lap_phi);
    conc = conc + (param.dtime)*(termc_final);
    
    
    if mod(istep,param.nprint)==0
        %display the iteration
        fprintf('Iteration = %d\n',istep);
        toc
        %display images
        subplot(2,2,1)
        imagesc(phi);
        title("phi")
        colorbar
        %colormap jet
        
        subplot(2,2,2)
        imagesc(theta)
        title("theta")
        colorbar
        %colormap jet
        
        subplot(2,2,3)
        imagesc(conc);
        title("conc")
        %caxis([0,1]);
        colorbar
        %colormap jet
        
        subplot(2,2,4)
        imagesc(tempr);
        title("tempr")
        colorbar
        %colormap jet
        drawnow
        
    end
    
end
end