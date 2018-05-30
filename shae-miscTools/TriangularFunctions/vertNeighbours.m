%%
function Vmat=vertNeighbours(TRI,neighbourhood)
% this function takes a triangulation and produces the neighbour matrix
% looking within a ring of size neighbourhood
if nargin==1
neighbourhood=1;
end
vertexNumber=max(TRI(:));
E=[TRI(:,[1,3]);TRI(:,2:3);TRI(:,1:2)];
[E,~,~]=unique(E,'rows');
%Eind=sub2ind([vertexNumber,vertexNumber],E(:,1),E(:,2));
Vmat=sparse(E(:,1),E(:,2),1,vertexNumber,vertexNumber,10*vertexNumber);
Vmat=(Vmat+Vmat');
Vmat=mpower2(Vmat,neighbourhood)>0;
% Vmat2 = Vmat;
% Vmat3 = Vmat;
% Vmat(diag(ones(size(Vmat,1),1))) = 0;
Vmat =Vmat.*(1-eye(size(Vmat)));

end
