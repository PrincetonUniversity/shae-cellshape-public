function [F2Vmat,invTri]=Face2Vert(TRI)
%Face2vert takes a triangulation and creates a matrix that takes the face
%indices to vertex indices. it also makes invTri, and matrix that tells
%which vertex indicies neighbour a given face.

vertexNumber=max(TRI(:));
faceNumber=length(TRI);
FTRI=1:length(TRI);
FTRI=repmat(FTRI',[1,3]);
pts=[FTRI(:),TRI(:)];


F2Vmat=sparse(pts(:,1),pts(:,2),1,faceNumber,vertexNumber);

if nargout==2
[row,col]=find(F2Vmat);

neighbournumber=sum(F2Vmat);

invTri=zeros(vertexNumber,max(neighbournumber));
for i=1:vertexNumber
    invTri(i,1:neighbournumber(i))=row(col==i);
end

end