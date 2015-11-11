%%%%%%%
%%%%%%%
% This code approximates the poles of the cells using the mean curvature.
% Each pole of the cell is considered to be the points on the surface of the
% center within a distance of 1/k from an endpoint of the cell, where k is
% the average mean curvature of the cell near that endpoint. 

%%%%%
% The poles are each represented using a mask, an array of points such that
% the value of the array at position (u,v) is one if (u,v) is in the pole
% and zero otherwise.


function [pole_left,pole_right]=pole(gm,slice_binary,CL)

distance_travelled = cumsum(sqrt(sum(diff(CL).^2,2)));

cell_length = distance_travelled(length(distance_travelled));

scale = 1.0;

radius_first = sum(nansum(slice_binary{1}))*scale/sum(nansum(slice_binary{1}.*gm));
radius_last = sum(nansum(slice_binary{length(CL)}))*scale/sum(nansum(slice_binary{length(CL)}.*gm));

firstpolecoordinate = find(distance_travelled > radius_first, 1, 'first');
lastpolecoordinate = find(cell_length-distance_travelled > radius_last, 1, 'last');

pole_left = zeros(size(slice_binary{1}));
pole_right = zeros(size(slice_binary{1}));

for i = 1:max(firstpolecoordinate,length(CL)-lastpolecoordinate);
    pole_left = pole_left+slice_binary{i};
    pole_right = pole_right + slice_binary{length(CL)-i+1};
end

end





