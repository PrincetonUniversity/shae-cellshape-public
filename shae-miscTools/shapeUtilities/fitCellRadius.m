
%%
function [r0,sigma,excessProb] = fitCellRadius(rPerGrid,areaPerGrid);
% This function will take in a surface and compute the "radius" of the
% cell. The model is that there is an object whose surface is paramaterized
% as distance from a centerline. There are three parameters to this model,
% the radius, relative magnitude of fluctuations and lengthscale over which
% the poles/endcaps exist. 

%%
 % for example, in triangles
% rPerGrid = Fcoord.r;
% [~,areaPerGrid,~]=triArea(F);
%%
% remove the sparse data points at the small side
lowCut = 20/nnz(areaPerGrid>0);
if nnz(areaPerGrid>0)<40
    r0 = nan;
    sigma = nan;
    excessProb = nan;
    return
end
minVal = quantile(rPerGrid(rPerGrid>0),lowCut);
rPerGrid(rPerGrid<=minVal) = nan;

% build a list of bins for histogramming
% use Freedman?Diaconis rule to determine bin widths
binWidth = 2*iqr(rPerGrid).*nnz(areaPerGrid>0)^(-1/3);
bins = minVal:binWidth:nanmax(rPerGrid(:));
% this can leave some radii without a home, so rebuild this many bins
% stretching from min to max.
nbins = length(bins);
nbins = max(nbins,100);
bins = linspace(minVal,nanmax(rPerGrid(:)),nbins);
[~,idx] = histc(rPerGrid,bins);

% those below the range, put at the top
idx(idx==0) = max(idx)+1;

% calculate probabilities and normalize
sum1 = accumarray(idx,areaPerGrid);
sum1 = sum1./sum(sum1);
% after accounting for those bins at the bottom, remove them
sum1(end) = [];
smallRadiusAreaFrac = 1-sum(sum1);

filtSigma = 5;
filterArray = exp(-((-25:25)).^2/(2*filtSigma.^2));
filterArray = filterArray./sum(filterArray);
filteredPdf = conv(sum1,filterArray,'same');
%% make a plot including raw and filtered data

%% start with some rough estimate of the peak width, then modify
lowEnd = round(0.05*nbins); % in bins
highEnd = round(0.05*nbins); % in bins

% find the maximum
[~,maxProbIdx] = max(filteredPdf);
if (maxProbIdx/nbins)<0.3;
    warning('cellRadius:LowPeak','Peak in probability near zero');
    maxProbIdx = 0.8*nbins;
end
% make sure the maximum is not too close to the end 
if (maxProbIdx+highEnd)>length(filteredPdf)
    highEnd = round(0.95*length(filteredPdf)-maxProbIdx);
end
% fit the region around the max to a parabola
bestFit1 = polyfit((maxProbIdx-lowEnd:maxProbIdx+highEnd)',...
    log(filteredPdf(maxProbIdx-lowEnd:maxProbIdx+highEnd)),2);
%%
% determine the "width" of the parabola in bins
Q = 1/bestFit1(1);

% from this first estimate, put some bounds on how wide for next round of
% fits
Q = min(abs(Q),nbins*0.4);
Q = max(abs(Q),4);
lowEnd = round(Q/4);
highEnd = round(Q/2);

% find the maximum
maxProbIdx = round(-bestFit1(2)/(2*bestFit1(1)));
% make sure the maximum is not too close to the end 
if (maxProbIdx+highEnd)>length(filteredPdf)
    maxProbIdx = length(filteredPdf)-highEnd-2;
end
% fit the region around the max to a parabola
bestFit1 = polyfit(bins(maxProbIdx-lowEnd:maxProbIdx+highEnd)',...
    log(filteredPdf(maxProbIdx-lowEnd:maxProbIdx+highEnd)),2);

% find the width of the parabola 
Q = 1/bestFit1(1);
sigma = sqrt(-Q/2);

% find the vertex of the parabola
r0 = -bestFit1(2)/(2*bestFit1(1));
y0 = bestFit1(3)-(r0^2/Q);

% evaluate this quadratic/guassian fit to determine the extra contribution 
% of small distances
peakEval = exp(polyval(bestFit1,bins));

% % % where do these cross and have non-negligible prob?)
% % crossingPt = find(and(...
% %     (filteredPdf-peakEval')<=0,...
% %     peakEval'>1e-10),1,'first');
% the difference in the observed pdf and the fit pdf is
% the fraction of surface area that is at a smaller radius than a cylinder
% would predict. We need to add to this the fraction of surface area for
% points that were at the smallest radius, where the sampling density is
% too low
if numel(filteredPdf)+1==numel(peakEval)
excessProb = sum(filteredPdf-transpose(peakEval(1:end-1)));
else
excessProb = sum(filteredPdf-peakEval');
end

if nargout == 0;
%% get ready to plot
figure();
% plot the raw data
plot(bins,sum1,'bx-');
hold on;
% plot the filtered data
plot(bins,filteredPdf,'r')
% plot the extrapolate fit
plot(bins,peakEval,'g:','linewidth',1);
% plot the fit in the range where it is being fit
plot(bins(maxProbIdx-lowEnd:maxProbIdx+highEnd),...
    exp(polyval(bestFit1,bins(maxProbIdx-lowEnd:maxProbIdx+highEnd))),...
    'g--','linewidth',2);
% add a title with the parameters of interest
title(['sigma = ',num2str(sigma),'.  x0 = ',num2str(r0),...
    '.  y0 = ',num2str(y0),'. Excess probability = ',...
    num2str(excessProb,'%4.3g'), ...
    ' + ', num2str(smallRadiusAreaFrac,'%4.3g'),... 
    ' = ', num2str(excessProb+smallRadiusAreaFrac,'%4.3g')]);
% set some reasonable display scales
ylim([1e-5,1e-1]);
xlimits = xlim;
set(gca,'yscale','log')
end

excessProb = excessProb+smallRadiusAreaFrac;
if r0<=0
    r0 = nan;
    excessProb = nan;
    sigma = nan;
end

if or(excessProb>1,excessProb<-0.1);
    r0 = nan;
    excessProb = nan;
    sigma = nan;
end
