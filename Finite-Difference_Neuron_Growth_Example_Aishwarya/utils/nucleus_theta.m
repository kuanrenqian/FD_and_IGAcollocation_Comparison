function [phi,tempr,theta,conct] = nucleus_theta(Nx,Ny,seed)

phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);
theta1 = rand(Nx,Ny);
conct = zeros(Nx,Ny);
% numx = 0.8*Nx;
% numy = 0.8*Ny;
% posx = randi(Nx,numx);
% posy = randi(Ny,numy);
% for i = 1:numx
%     for j = 1:numy
%         theta1(posx(i),posy(j))=rand(1,1);
%     end
%end

%conc = 0.75.*ones(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        incl = atan2((i-Ny/2),(j-Nx/2));
        if(i>=Nx/2 && j>=Ny/2)
            %incl = -incl;
            theta1(i,j) = 0.5+ 0.5*incl/(pi);  
        end
        
        if ((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed)
            phi(i,j) = 1.0;
            r = sqrt((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2));
            phi(i,j) = 1.0;
            conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
        end
    end
end

theta = theta1;