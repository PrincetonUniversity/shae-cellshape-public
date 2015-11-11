function [Area,Amesh,V]=surfArea(x,y,z)

%Takes a closed contour of the form S=(x(u,v),y(u,v),z(u,v)) and calculates
%the surface area. 
%outputs Area, Area mesh for each pixel, and normalized area vector V

  [Xu(:,:,1),Xv(:,:,1)]=gradient(x);
    [Xu(:,:,2),Xv(:,:,2)]=gradient(y);
    [Xu(:,:,3),Xv(:,:,3)]=gradient(z);

            A=(cross(Xu,Xv,3));
             Amesh=sqrt(sum(A.^2,3));
            Area=sum(sum(Amesh));           
A(:,1,:)=repmat(mean(A(:,2,:)),[size(A,1),1,1]);
A(:,end,:)=repmat(mean(A(:,end-1,:)),[size(A,1),1,1]);

Amag=sqrt(sum(A.^2,3));
V=bsxfun(@rdivide,A,Amag);

%fix orientations so that normals point on average, away form the center of
%mass
V1=reshape(V,[],3);
x=x(:)-mean(x(:));
y=y(:)-mean(y(:));
z=z(:)-mean(z(:));
%%
cmV=[x,y,z];
cmV=normr(cmV);
if sum(dot(cmV,V1,2))<0
    V=-V;
end


            
