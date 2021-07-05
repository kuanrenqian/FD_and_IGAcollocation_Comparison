function [epsilon, epsilon_deriv, aap, P_dy, P_dx] = kqGetEpsilonAndAap(epsilonb,delta,phi,theta,NuNv,NuN1v,N1uNv)
% kqGet_epsilon_and_aap(phi,theta,NuN1v,N1uNv) 
% This function calculates epsilon and aap (a*a') based on phi, theta,
% NuN1v, and N1uNv.
% kqGet_epsilon_and_aap() output epsilon aap in the format of [epsilon,aap]

aniso = 6;

P_dy = NuN1v*phi;
P_dx = N1uNv*phi;

sz = sqrt(length(P_dx));
% P_dy = reshape(P_dy,sz,sz);
% P_dx = reshape(P_dx,sz,sz);
% for j = 1:sz
%     for i = 1:floor(sz/2)
%         P_dx(j,i) = abs(P_dx(j,i));
%         P_dy(i,j) = abs(P_dy(i,j));
%     end
%     for i = floor(sz/2)+1:sz
%         P_dx(j,i) = -abs(P_dx(j,i));
%         P_dy(i,j) = -abs(P_dy(i,j));
%     end
% end      

P_dy = reshape(P_dy,sz*sz,1);
P_dx = reshape(P_dx,sz*sz,1);

% atheta =smooth(full(atan2(phidy,phydx)));
atheta =(full(atan2(P_dy,P_dx)));
xtheta = 2*pi*(NuNv*theta);
% xtheta = 0.2;

epsilon = epsilonb.*(1.0+delta*cos(aniso*(atheta-xtheta)));
epsilon_deriv = -epsilonb.*(aniso*delta*sin(aniso.*(atheta-xtheta)));
aap = epsilon.*epsilon_deriv;

epsilon = NuNv\epsilon;
aap = NuNv\aap;

epsilon = sparse(epsilon);
epsilon_deriv = sparse(epsilon_deriv);
aap = sparse(aap);
