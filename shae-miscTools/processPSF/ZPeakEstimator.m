function [bestEst,output] = bpbZPeakEstimator(pkCutOut,param,zPos,options)

% % param is the same structure as overall params

% % options have specific parameters/options
% options.intensityMode -> what to use as 'intensity'

% choose what 'intensity' function to use
switch options.intensityMode
    case 'Brenner'
        intensity = squeeze(sum(sum(...
            (pkCutOut(1:end-2,:,:)-pkCutOut(3:end,:,:)).^2 ...
            ,1),2));
    case 'Sum'
        intensity = squeeze(sum(sum(pkCutOut,1),2));
        
end

[~,maxP] = max(intensity);
fitRadius = round(abs(param.zRadius/((zPos(2)-zPos(1)))));

minZ = max(1,maxP-fitRadius);
maxZ = min(length(zPos),maxP+fitRadius);

fitRegionX = zPos(minZ:maxZ);
intensityFit = intensity(minZ:maxZ);

[aOut] = polyfit(fitRegionX',intensityFit,param.nPoly);

zerosA = roots(aOut(1:param.nPoly).*(param.nPoly:-1:1));
zerosA(100*abs(imag(zerosA))>abs(real(zerosA))) = [];

zerosA(((real(zerosA)))>max(zPos([maxZ,minZ]))) = [];
zerosA(((real(zerosA)))<min(zPos([maxZ,minZ]))) = [];
zerosVal = polyval(aOut,zerosA);
[~,maxId] = max(real(zerosVal));
bestEst = (zerosA(maxId));

maxInt = polyval(aOut,(bestEst));


% find local peaks
isAboveA = 1;
isAboveB = 1;
countsA = 0;
countsB = 0;
for iiZ = maxP:size(intensity,1);
    
    % less than 15% to above 25%
    if isAboveA
        if intensity(iiZ)<param.tempA1*maxInt
            isAboveA= 0;
        end
    else
        if intensity(iiZ)>param.tempA2*maxInt
            isAboveA = 1;
            countsA = countsA+1;
        end
    end
    
    % less than 25% to above 50%
    if isAboveB
        if intensity(iiZ)<param.tempB1*maxInt
            isAboveB= 0;
        end
    else
        if intensity(iiZ)>param.tempB2*maxInt
            isAboveB = 1;
            countsB = countsB+1;
        end
    end
end

isAboveA = 1;
isAboveB = 1;
for iiZ = maxP:-1:1
    % less than A1 to above A2
    if isAboveA
        if intensity(iiZ)<param.tempA1*maxInt
            isAboveA= 0;
        end
    else
        if intensity(iiZ)>param.tempA2*maxInt
            isAboveA = 1;
            countsA = countsA+1;
        end
    end
    
    % less than B1 to above B2
    if isAboveB
        if intensity(iiZ)<param.tempB1*maxInt
            isAboveB= 0;
        end
    else
        if intensity(iiZ)>param.tempB2*maxInt
            isAboveB = 1;
            countsB = countsB+1;
        end
    end
end

%  display plots
if param.displayIntensity
    %%
    figure(gcf);
    clf;
    plot(zPos,intensity,'bo:');
    hold on;
    plotFitX = linspace(fitRegionX(1),fitRegionX(end),100);
    plot(plotFitX ,polyval(aOut,(plotFitX)),'k');
    ylimits = ylim;
    plot([zPos(end),zPos(1)],[param.tempA1.*maxInt,param.tempA1.*maxInt],'k:');
    plot([zPos(end),zPos(1)],[param.tempA2.*maxInt,param.tempA2.*maxInt],'k:');
    plot([zPos(end),zPos(1)],[param.tempB1.*maxInt,param.tempB1.*maxInt],'k:');
    plot([zPos(end),zPos(1)],[param.tempB2.*maxInt,param.tempB2.*maxInt],'k:');
    
    ylabel(options.intensityMode);
    axis tight;
    xlabel('stage position');
end

% if all possible peaks in this zStack have been rejected, ignore this
% zStack
if isempty(bestEst)
    bestEst = -inf;
    
end

isRejected = or(or(... % reject
    (countsA+countsB)>param.tempA,...
    (bestEst<(min( zPos(round(fitRadius/2)),...
    zPos(end-round(fitRadius/2)))))),...
    max(intensity)<param.zIntensityMin...
    );



if isRejected
    if param.displayIntensity
        plot([bestEst,bestEst],ylimits,'r');
    end
    bestEstP = nan;
    fwhmZ = nan;
    fwhmZ1 = nan;
    fwhmZ2 = nan;
    
else % keep and calculate other properties
    if param.displayIntensity
        plot([bestEst,bestEst],ylimits,'g');
    end
    
    % % full width at half max in Z along central portion
    belowHalfInt = intensity<0.5*maxInt;
    % first estimate is to the pixel level
    fwhmZ1 = find(belowHalfInt(1:maxP),1,'last');
    fwhmZ2 = find(belowHalfInt(maxP:end),1,'first')+maxP-1;
    % second estimate is a linear approximation between neighboring pixels
    try
        fwhmZ1A = ((0.5*maxInt-intensity(fwhmZ1)).*...
            (zPos(fwhmZ1-1)-zPos(fwhmZ1))./(intensity(fwhmZ1-1)-intensity(fwhmZ1)))+...
            zPos(fwhmZ1);
        fwhmZ2A = ((0.5*maxInt-intensity(fwhmZ2)).*...
            (zPos(fwhmZ2+1)-zPos(fwhmZ2))./(intensity(fwhmZ2+1)-intensity(fwhmZ2)))+...
            zPos(fwhmZ2);
        fwhmZ = abs(fwhmZ2A-fwhmZ1A);
        
        
        
    catch ME
        fwhmZ = nan; % runs into an error in fwhm determination
        isRejected = 1;
    end
    
end

if param.displayIntensity
    if isempty(param.dispPause);
        
        pause();
    else
        pause(param.dispPause);
    end
end

if isempty(fwhmZ)
    fwhmZ = nan;
    isRejected = 1;
end



% verbose output
if nargout>1
    % intensity type specific
    output.aOut = aOut;
    output.intensity = intensity;
    output.maxInt = maxInt;
    output.bestEst = bestEst;
    output.isKept = not(isRejected);
    output.fwhmZ = fwhmZ;
    output.fwhmZ1 = fwhmZ1;
    output.fwhmZ2 = fwhmZ2;
    
    % error
    if exist('ME');
        output.ME = ME;
    end
end