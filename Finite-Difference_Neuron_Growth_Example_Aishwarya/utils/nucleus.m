function [phi,tempr,conct] = nucleus_old(Nx,Ny,seed)

phi = zeros(Nx,Ny);
tempr = zeros(Nx,Ny);
conct = zeros(Nx,Ny);
for i=1:Nx
    for j=1:Ny
        if ((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed)
            r = sqrt((i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2));
            phi(i,j) = 1.0;
            conct(i,j) = (0.5+0.5*tanh((sqrt(seed)-r)/2));
        end
    end
end

