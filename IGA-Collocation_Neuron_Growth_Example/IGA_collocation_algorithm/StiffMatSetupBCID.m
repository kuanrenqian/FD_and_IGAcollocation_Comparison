function [coll_Lhs, coll_Rhs] = StiffMatSetupBCID(coll_Lhs, coll_Rhs,bcid,N)

bcdof = find(bcid==1);
diag_lhs = spdiags(coll_Lhs,0).*(1-bcid)+bcid;
coll_Lhs(bcdof,:) = 0;
coll_Lhs(:,bcdof) = 0;
coll_Rhs(bcdof) = N(bcdof);
coll_Lhs = spdiags(diag_lhs, 0, coll_Lhs);

coll_Lhs = sparse(coll_Lhs);
coll_Rhs = sparse(coll_Rhs);