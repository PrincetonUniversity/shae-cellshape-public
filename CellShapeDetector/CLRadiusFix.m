
function [radius,CLinput,slice_binary,avg_radius]=CLRadiusFix(xf,yf,zf,CL)
% % BPB notes 20130911
% - input x, y, z are parameterized as u hoops, each with v elements
% - the 'radius' for each surface element is measured to the closest point 
% on an interpolated centerline, not necessarily one that has the same u 
% coordinate
% - the average hoop radius is then the average radius for a given u,
% weighted by the surface area of that element
% - the output centerline is the same as the input centerline. There was
% some attempt to make each centerline point U minimize the distance to the
% surface points in that U hoop, but this results in someone non-uniformly
% spaced points along the centerline.

% number of points added between points in CL in the interpolation
% %  BPB side comment: this is not quite correct, there will be N times as many points in the
% % interpolated centerline and the first and last points will be the same.
% % The intermediate points are not necessarily in the interpolated contour
% % eg: 
% % Start with 3 points, N=4, final number of points, m = 12
% % Start with 3 points, insert four between each, m = 11 = 3 (intial)+4*(two segments)
N = 10; 

% calculate the size of U and V
sizes = size(xf);

% surface area for weighting
[~, Amesh, ~] = surfArea(xf,yf,zf);

% temp holding place for the centerline coordinates, return these
% unperturbed as the output centerline.
CLinput = CL;

% make surface points into an array
SurfPointsvect = cat(2, xf(:), yf(:), zf(:));

% U indices to interpolate centerline
interpCLpoints = linspace(1,size(CL,1),N*size(CL,1));

% interpolated centerline
CLprime = cat(2,interp1(double(CL(:,1)),interpCLpoints)',...
    interp1(double(CL(:,2)),interpCLpoints)',...
    interp1(double(CL(:,3)),interpCLpoints)');

% find the shortest distance between the surface points and the centerline
[distances,indices] = pdist2(CLprime,SurfPointsvect,'minkowski',2,'Smallest',1); 

% reshape this into the UV shape
radius = double(reshape(distances, sizes(1), sizes(2)));


avg_radius = sum(radius.*Amesh)./sum(Amesh);

% % % for each column or U coordinate, find the centerline point that minimizes
% % % the distances to the hoop assigned
% % for iU = 1:sizes(2)
% %     % fragment of the centerline, interpolated between iU and iU+1
% % %     tempCL = CLprime((1:N)+(N*(iU-1)),:);
% %     % find all the pairwise distances
% %     distances = pdist2(CLprime,[xf(:,iU),yf(:,iU),zf(:,iU)],'euclidean');
% %     % weight by appropriate surface areas of section
% %     distances = bsxfun(@times, Amesh(:,iU)',distances);
% %     % find the smallest sum
% %     [sumDist, location] = min(sum(distances,2));
% %     % convert to average
% %     avg_radius(iU) = sumDist/(sum(Amesh(:,iU)));
% %     % change the position of the Uth centerline point to this one
% %     CL(iU,:) = CLprime(location,:);
% % end
% % the final point and initial point sometimes get confused about where to
% % be because of the nans and singularities
% % CL(end,:) = CLinput(end,:);
% % CL(1,:) = CLinput(1,:);

% the slice_binary is used in "pole" something about determining poles
 indices1 = double(reshape(indices, sizes(1), sizes(2)));
 numberslices = length(CL); 
 slice_binary = cell(1, numberslices);
 divisor = ceil((max(indices)- min(indices))/numberslices);
 maxcount = 0;
% 
% %%% Partitions the surface into slices
for i=1:1:(numberslices)
    %
    mincount = maxcount;
    maxcount = mincount + divisor;
    %
    slice_binary{i} = +((indices1 > mincount) & (indices1 <= maxcount));
end
% slice_binary = nan;