%% Run FitDataCompareTri.m on all the folders
function [] = enrichmentSmoothingSpline();

% clear all
FolderList = uipickfiles('Prompt','Select folders you want for curvature analysis');
% N3 = enrichment numbers
% counts = bottom histogram
% X = curvature value (x/um2)
% N3se = error
%%
% these come from FitDataCompare code
% control=varargin{1};
% controlf=varargin{2};
% path=varargin{3};
% options=varargin{4};
options.undersample=0; % 1 for normal, 0 for super resolution
options.midflag=0; % leave on 0 to use all the 3D data
options.poleflag=1;% 1 leaves poles in, 0 throws out poles after normal. -1 before norm
options.plotflag=0;
% options.backgroundSubtract = -999; %use min to max of data
% options.backgroundSubtract = 150;
options.backgroundSubtract = -999;


progressbar('tic toc');
for iFolder = 1:length(FolderList);
    if options.plotflag
        figure()
    end
    [N3,counts,X,N3se,g,gg,NN,N2] = FitDataCompareTri('gc','Ff',FolderList{iFolder},options);
    
    progressbar(iFolder/length(FolderList));
end

%% Start here given enrichment data

FileList = uipickfiles('Filterspec','*curve.mat');

%%

% exclude data threshold (curvatures with low representation), don't use bins without many samples
threshold = 1e-3;

fig1 = figure(); % Enrichment splines and error estimates
fig2 = figure(); % histograms of available curvature
fig3 = figure(); % fraction of total signal

% set limits for where to display the data
xlimits = [-20,30];

% smoothing class, how hard to smooth the data for curvature probabilities
smoothingClass = 'very heavy'; % 'none','default','light','moderate','heavy','very heavy'


% choose a list of colors
colorsList = parula(length(FileList));

for iFile = 1:length(FileList);
    % load in the files one at a time
    tempData = load(FileList{iFile});
    
    % attempt to choose a scale for smoothing based on the data, whether it
    % is Gaussian or K1, or centerline curvature, etc, different magnitudes
    if abs(tempData.Fittingdata.X(2)-tempData.Fittingdata.X(3))<1e-2
        tempData.Fittingdata.X = tempData.Fittingdata.X*1e3;
    end
    
    % smooth data using a smoothing spline
    [fitresult, gof] = createFit(tempData.Fittingdata.X, tempData.Fittingdata.N3,tempData.Fittingdata.counts, threshold,0.97);
    
    % postpend with nans to get to the correct length
    while length(tempData.Fittingdata.N3se)<length(tempData.Fittingdata.N3)
        tempData.Fittingdata.N3se(end+1) = nan;
    end
    
    % smooth error estimates with a smoothing spline
    [errFit, gofErr] = createFit(tempData.Fittingdata.X, tempData.Fittingdata.N3se,tempData.Fittingdata.counts, threshold);
    
    % create a temporary figure to grab the xData
    figTemp = figure();
    h = plot(fitresult);
    xVals = get(h,'xData');
    close(figTemp);
    
    % plot the result on the main figure
    figure(fig1);
    h = plot(xVals,fitresult(xVals));
    fileSpecificColor = colorsList(iFile,:);
    set(h,'DisplayName',tempData.Fittingdata.label, 'color', fileSpecificColor, 'linewidth',2);
    
    % Label axes
    xlabel( 'Geometric parameter' );
    ylabel( 'Enrichment' );
    hold on
    
    % display the error estimates
    h3(1) = plot(xVals,fitresult(xVals)+errFit(xVals));
    h3(2) = plot(xVals,fitresult(xVals)-errFit(xVals));
    
    set(h3,'DisplayName',[tempData.Fittingdata.label,'_error'],'color',fileSpecificColor);
    
    
    % plot the histograms of curvature values
    figure(fig2);
    subplot(length(FileList),1,iFile);
    bar(tempData.Fittingdata.X,tempData.Fittingdata.counts,...
        'EdgeColor',fileSpecificColor,'FaceColor','none');
    
    % csaps smoothing occurs when the smoothing parameter p is set based on
    % the spacing of data
    h = mean(diff(tempData.Fittingdata.X));
    switch smoothingClass
        case 'none'
            p = 1;
        case 'light'
            p = 1/(1 + (h^3)/60);
        case {'default','moderate'}
            p = 1/(1 + (h^3)/6);
        case 'heavy'
            p = 1/(1 + (h^3)/0.6);
        case 'very heavy'
            p = 1/(1 + (h^3)/0.06);
        otherwise  % default
            p = 1/(1 + (h^3)/6);
    end
    % calculate smoothing spline
    pCurvatureSpline = csaps(tempData.Fittingdata.X,tempData.Fittingdata.counts,p);
    pCurvatureXVals = linspace(tempData.Fittingdata.X(1),tempData.Fittingdata.X(end),1000);
    
    % display smoothing spline
    hold on;
    plot(pCurvatureXVals,fnval(pCurvatureSpline,pCurvatureXVals),':','color',fileSpecificColor,'linewidth',2);
    
    xlim(xlimits);
    
    ylabel('fraction of surface area');
    xlabel('gaussian curvature (\mum^2)');
    
    xVals1 = tempData.Fittingdata.X;
    xValMin = find(tempData.Fittingdata.X<xVals(1),1,'last');
    xValMax = find(tempData.Fittingdata.X>xVals(end),1,'first');
    
    % plot the fraction of signal from each bin
    figure(fig3);
    subplot(length(FileList),1,iFile);
    
    % calculate probability and intensity
    aggregatedResultsInWorkspace(iFile).X = tempData.Fittingdata.X(xValMin:xValMax);
    aggregatedResultsInWorkspace(iFile).P = tempData.Fittingdata.counts(xValMin:xValMax);
    aggregatedResultsInWorkspace(iFile).IoverIstar = tempData.Fittingdata.N3(xValMin:xValMax);
    
    % display the result
    plot(tempData.Fittingdata.X(xValMin:xValMax),...
        fitresult(tempData.Fittingdata.N3(xValMin:xValMax)).*tempData.Fittingdata.counts(xValMin:xValMax),...
        'Color',fileSpecificColor);
    xlim(xlimits);
    
    % label axes
    title(num2str(nansum(...
        fitresult(tempData.Fittingdata.N3(xValMin:xValMax)).*tempData.Fittingdata.counts(xValMin:xValMax))));
    ylabel('average value * fraction of surface area');
    xlabel('gaussian curvature (\mum^2)');
end

% switch back to the curvature enrichment figure
figure(fig1);
h2 = legend('Location', 'best');
set(h2,'interpreter','none' );

xlim(xlimits);
ylim([0.1,1.3]);
plot(xlim,[1,1],'color', 'k', 'linewidth',4,'DisplayName','uniform coverage')
%%


%%




