function [Area,vertArea,faceArea]=triArea(F)
% function takes patch object F and calculates total area of object, area
% at each vertex, and area at each face. 

fxyz=F.vertices(F.faces(:),:);
fxyz=reshape(fxyz,[],3,3);

V1=squeeze(fxyz(:,1,:)-fxyz(:,2,:));
V2=squeeze(fxyz(:,1,:)-fxyz(:,3,:));
V3=squeeze(fxyz(:,2,:)-fxyz(:,3,:));
faceArea=cross(V1,V2,2)/2;
faceArea=sqrt(sum(faceArea.^2,2));
vertArea=accumarray(F.faces(:),[faceArea;faceArea;faceArea],[],@sum)/3;
Area=sum(faceArea);