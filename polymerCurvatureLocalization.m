function [ rcCurvature ] = polymerCurvatureLocalization( curvature,rc )
%POLYMERCURVATURELOCALIZATION takes a polymer defined from the polyFit
% structure and interpolates to find the average curvature at those points
% rc(:,1)=mod(rc(:,1)-1,size(curvature,1));
% rc(:,2)=mod(rc(:,2)-1,size(curvature,2));
% rc=rc+1;
curvature=[curvature;curvature;curvature];
rcCurvature=(interp2(curvature,rc(:,1),rc(:,2)));
end

