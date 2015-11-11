function outputColors = jetkw(nColors);
% this function takes the original jet colorscale and makes the max color
% white and the min color black
if nargin<1
   nColors = size(get(gcf,'colormap'),1);

end
tempColor = jet(nColors);

%graymap lookup
grayscaleIndex = linspace(0,1,nColors)-0.5;
grayscaleIndex = grayscaleIndex.^(7);
grayscaleIndex = grayscaleIndex(:);
maxVal = max(abs(grayscaleIndex(:)));
grayscaleIndex = grayscaleIndex./maxVal;
grayscaleOffset = ones(nColors,3);
grayscaleOffset = bsxfun(@times,grayscaleOffset,grayscaleIndex);
% augment by grays
tempColor = tempColor+grayscaleOffset;
% make sure to stay in valid regime
tempColor = max(tempColor,0);
outputColors = min(tempColor,1);



