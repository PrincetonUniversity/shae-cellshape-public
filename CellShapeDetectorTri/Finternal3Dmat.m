function [FL,FU]=Finternal3Dmat(F,cline_para,option)
if nargin==2
    option='uniform';
end

% create upper and lower triangular matricies for the internal energy of a
% rectangular mesh of size x, using parameters in clinepara. 

Fmat=vertNeighbours(F.faces,1);
Fmat=bsxfun(@rdivide,Fmat',sum(Fmat,2))';


    K=size(Fmat,2);
nMatrix=Fmat-speye(size(Fmat));
nMatrix(isnan(nMatrix))=0;
nMatrix=nMatrix';
f=(cline_para.gamma*speye(K)-cline_para.alpha3d*nMatrix+cline_para.beta3d*nMatrix*nMatrix);
[FL,FU]=lu(f);
%force=inv(cline_para.gamma*speye(K,K)-cline_para.alpha*nMatrix+cline_para.beta*nMatrix*nMatrix);

