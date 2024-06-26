function [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,size_collpts] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs)
% kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv)
% This function calculates all basis component for IGA-collocation 2D structured grid
% kqCollocationDers() takes 2 knotvectors in u and v direction with corresponding
% curve order for calculation, order_deriv means the order of derivative
% that you want to calculate.
% Note that kqCollocationDers() automatically scales input knotvector using
% max(knotvector)

lenu = length(knotvectorU)-2*(p-1);
lenv = length(knotvectorV)-2*(q-1);

% scale knotvector to 0~1 - NURBS book convention
% knotvectorU = knotvectorU/max(knotvectorU);
% knotvectorV = knotvectorV/max(knotvectorV);

k=1;
coll_p = zeros(length(knotvectorU)-p,2);
for i = 1:length(knotvectorU)-p-1
    coordx = sum(knotvectorU(i+1:i+p))/p;
    for j = 1:length(knotvectorV)-q-1
        coordy = sum(knotvectorV(j+1:j+q))/q;
        coll_p(k,:) = [coordx,coordy];
        k=k+1;
    end
end
size_collpts = sqrt(length(coll_p));

ind_i = [];
ind_j = [];

NuNv_val = [];
N1uNv_val = [];
NuN1v_val = [];
N1uN1v_val = [];
N2uNv_val = [];
NuN2v_val = [];
N2uN2v_val = [];

% number of basis function (setup this way*) (current version works)
% using a modified version of derbasisfun3, need to fix this!!!!!
nobu = length(knotvectorU)-p-2; 
nobv = length(knotvectorV)-p-2;

knotSpanU = zeros([size_collpts^2,1]);
dersU = zeros([3,p+1,size_collpts^2]);
knotSpanV = zeros([size_collpts^2,1]);
dersV = zeros([3,q+1,size_collpts^2]);
ksU = zeros([size_collpts^2,p+1]);
ksV = zeros([size_collpts^2,q+1]);

for i = 1:size_collpts
    for j = 1:size_collpts
        k = (i-1)*size_collpts+j;
        % extract u,v position
        u = coll_p(k,1);
        v = coll_p(k,2);

        % calculating basis function, its 1st & 2nd derivatives based on
        % collocation points and knot vector
        knotSpanU(k) = FindSpan_modified(nobu,p,u,knotvectorU);
        dersU(:,:,k)= dersbasisfuns3_modified(knotSpanU(k),p,order_deriv,u,knotvectorU);
        knotSpanV(k) = FindSpan_modified(nobv,q,v,knotvectorV);
        dersV(:,:,k)= dersbasisfuns3_modified(knotSpanV(k),q,order_deriv,v,knotvectorV);

        % extracting basis value and second derivatives from dersU&V
        Nu = dersU(1,:,k);
        Nv = dersV(1,:,k);
        N1u = dersU(2,:,k);
        N1v = dersV(2,:,k);
        N2u = dersU(3,:,k);
        N2v = dersV(3,:,k);

        % knot span vector -q/p -> for locating corresponding T
        ksU(k,:) = (knotSpanU(k):(knotSpanU(k)+p))-p;
        ksV(k,:) = (knotSpanV(k):(knotSpanV(k)+q))-q;

        for l = 1:length(ksU(k,:))
            for m = 1:length(ksV(k,:))
                ind = ksU(k,l)*lenv+ksV(k,m)+1; % calculating corresponding position using knotspan
                ind_i(end+1) = k;
                ind_j(end+1) = ind;
                NuNv_val(end+1) = Nu(l)*Nv(m); % assigning Ni*Nj value
                N1uNv_val(end+1) = N1u(l)*Nv(m); % assigning N1i*Nj value
                NuN1v_val(end+1) = Nu(l)*N1v(m); % assigning Ni*N1j value
                N1uN1v_val(end+1) = N1u(l)*N1v(m); % assigning N1i*N1j value
                N2uNv_val(end+1) = N2u(l)*Nv(m); % assigning N2i*Nj value
                NuN2v_val(end+1) = Nu(l)*N2v(m); % assigning Ni*N2j value
                N2uN2v_val(end+1) = N2u(l)*N2v(m); % assigning N2i*N2j value
            end
        end
    end
end

NuNv = sparse(ind_i,ind_j,NuNv_val,length(coll_p),lenu*lenv);
N1uNv = sparse(ind_i,ind_j,N1uNv_val,length(coll_p),lenu*lenv);
NuN1v = sparse(ind_i,ind_j,NuN1v_val,length(coll_p),lenu*lenv);
N1uN1v = sparse(ind_i,ind_j,N1uN1v_val,length(coll_p),lenu*lenv);
N2uNv = sparse(ind_i,ind_j,N2uNv_val,length(coll_p),lenu*lenv);
NuN2v = sparse(ind_i,ind_j,NuN2v_val,length(coll_p),lenu*lenv);
N2uN2v = sparse(ind_i,ind_j,N2uN2v_val,length(coll_p),lenu*lenv);

%%

% function [NuNv,N1uNv,NuN1v,N1uN1v,N2uNv,NuN2v,N2uN2v,coll_p,size_collpts,Control_points,dersU] = kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv,sprs)
% % kqCollocationDers(knotvectorU,p,knotvectorV,q,order_deriv)
% % This function calculates all basis component for IGA-collocation 2D structured grid
% % kqCollocationDers() takes 2 knotvectors in u and v direction with corresponding
% % curve order for calculation, order_deriv means the order of derivative
% % that you want to calculate.
% % Note that kqCollocationDers() automatically scales input knotvector using
% % max(knotvector)
% 
% lenu = length(knotvectorU)-2*(p-1);
% lenv = length(knotvectorV)-2*(q-1);
% 
% % scale knotvector to 0~1 - NURBS book convention
% % knotvectorU = knotvectorU/max(knotvectorU);
% % knotvectorV = knotvectorV/max(knotvectorV);
% 
% k=1;
% coll_p = zeros(length(knotvectorU)-p,2);
% for i = 1:length(knotvectorU)-p-1
%     coordx = sum(knotvectorU(i+1:i+p))/p;
%     for j = 1:length(knotvectorV)-q-1
%         coordy = sum(knotvectorV(j+1:j+q))/q;
%         coll_p(k,:) = [coordx,coordy];
%         k=k+1;
%     end
% end
% size_collpts = sqrt(length(coll_p));
% 
% %Control points
% size_cp = size_collpts-2;
% Control_points = zeros(size_cp^2,2);
% k=1;
% for i = 1:size_cp
%     for j = 1:size_cp
%         Control_points(k,1) = knotvectorU(i+p);
%         Control_points(k,2) = knotvectorV(j+q);
%         k = k+1;
%     end
% end
% 
% NuNv = zeros([length(coll_p),lenu*lenv]);
% N1uNv = NuNv;
% NuN1v = NuNv;
% N1uN1v = NuNv;
% N2uNv = NuNv;
% NuN2v = NuNv;
% N2uN2v = NuNv;
% 
% % number of basis function (setup this way*) (current version works)
% % using a modified version of derbasisfun3, need to fix this!!!!!
% nobu = length(knotvectorU)-p-2; 
% nobv = length(knotvectorV)-p-2;
% 
% knotSpanU = zeros([size_collpts^2,1]);
% dersU = zeros([3,p+1,size_collpts^2]);
% knotSpanV = zeros([size_collpts^2,1]);
% dersV = zeros([3,q+1,size_collpts^2]);
% ksU = zeros([size_collpts^2,p+1]);
% ksV = zeros([size_collpts^2,q+1]);
% 
% for i = 1:size_collpts
%     for j = 1:size_collpts
%         k = (i-1)*size_collpts+j;
%         % extract u,v position
%         u = coll_p(k,1);
%         v = coll_p(k,2);
% 
%         % calculating basis function, its 1st & 2nd derivatives based on
%         % collocation points and knot vector
%         knotSpanU(k) = FindSpan_modified(nobu,p,u,knotvectorU);
%         dersU(:,:,k)= dersbasisfuns3_modified(knotSpanU(k),p,order_deriv,u,knotvectorU);
%         knotSpanV(k) = FindSpan_modified(nobv,q,v,knotvectorV);
%         dersV(:,:,k)= dersbasisfuns3_modified(knotSpanV(k),q,order_deriv,v,knotvectorV);
% 
%         % extracting basis value and second derivatives from dersU&V
%         Nu = dersU(1,:,k);
%         Nv = dersV(1,:,k);
%         N1u = dersU(2,:,k);
%         N1v = dersV(2,:,k);
%         N2u = dersU(3,:,k);
%         N2v = dersV(3,:,k);
% 
%         % knot span vector -q/p -> for locating corresponding T
%         ksU(k,:) = (knotSpanU(k):(knotSpanU(k)+p))-p;
%         ksV(k,:) = (knotSpanV(k):(knotSpanV(k)+q))-q;
% 
%         for l = 1:length(ksU(k,:))
%             for m = 1:length(ksV(k,:))
%                 ind = ksU(k,l)*lenv+ksV(k,m)+1; % calculating corresponding position using knotspan
%                 NuNv(k,ind) = Nu(l)*Nv(m); % assigning Ni*Nj value
%                 N1uNv(k,ind) = N1u(l)*Nv(m); % assigning N1i*Nj value
%                 NuN1v(k,ind) = Nu(l)*N1v(m); % assigning Ni*N1j value
%                 N1uN1v(k,ind) = N1u(l)*N1v(m); % assigning N1i*N1j value
%                 N2uNv(k,ind) = N2u(l)*Nv(m); % assigning N2i*Nj value
%                 NuN2v(k,ind) = Nu(l)*N2v(m); % assigning Ni*N2j value
%                 N2uN2v(k,ind) = N2u(l)*N2v(m); % assigning N2i*N2j value
%             end
%         end
% 
% %         sNuNv = sum(NuNv(k,:));
% %         sN1uNv = sum(N1uNv(k,:));
% %         sNuN1v = sum(NuN1v(k,:));
% %         sN2uNv = sum(N2uNv(k,:));
% %         sNuN2v = sum(NuN2v(k,:));
% %         sN2uN2v = sum(N2uN2v(k,:));
% % 
% %         NuNv(k,:) = NuNv(k,:)/sNuNv;
% %         N1uNv(k,:) = N1uNv(k,:)/sNuNv - NuNv(k,:)*sN1uNv/(sNuNv^2);
% %         NuN1v(k,:) = NuN1v(k,:)/sNuNv - NuNv(k,:)*sNuN1v/(sNuNv^2);
% % 
% %         N2uNv(k,:) = N2uNv(k,:)/sNuNv - (2*N1uNv(k,:)*sN1uNv+NuNv(k,:)*sN2uNv)/(sNuNv^2) ...
% %             + (2*NuNv(k,:)*sN1uNv^2)/(sNuNv^3);
% %          NuN2v(k,:) = NuN2v(k,:)/sNuNv - (2*NuN1v(k,:)*sNuN1v+NuNv(k,:)*sNuN2v)/(sNuNv^2) ...
% %             + (2*NuNv(k,:)*sNuN1v^2)/(sNuNv^3);       
% % 
% %         N2uN2v(k,:) = N2uNv(k,:)/sNuNv - (N1uNv(k,:)*sNuN1v+NuN1v(k,:)*sN1uNv+NuNv(k,:)*sN2uN2v)/(sNuNv^2) ... 
% %             +(2*NuNv(k,:)*sN1uNv*sNuN1v)/(sNuNv^3);
%     end
% end
% 
% if(sprs ==1)
%     NuNv = sparse(NuNv);
%     N1uNv = sparse(N1uNv);
%     NuN1v = sparse(NuN1v);
%     N1uN1v = sparse(N1uN1v);
%     N2uNv = sparse(N2uNv);
%     NuN2v = sparse(NuN2v);
%     N2uN2v = sparse(N2uN2v);
% end