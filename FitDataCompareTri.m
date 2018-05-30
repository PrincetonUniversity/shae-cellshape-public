function [N3,counts,X,N3se,g,gg,NN,N4]= FitDataCompareTri(varargin)

%%%%%%%%%
%[N3,counts,X,N3se,g,gg]= FitDataCompare(varargin) creates a plot of a cell shape
%metric and the corresponding fluorescent signal in the bins of that
%histogram. N3 is the relative enrichment, counts are the counts of each
%curvature bin, X is the curvature bins, N3se is the standard error on the
%erichment signal, g and gg are the total shape and signal values.
% Inputs:
%    FitDataCompare -  select folder of fit files and compare gaussian
%    curvature with fluorescent signal.
%    FitDataCompare(ShapeMetric) - compare some shapemetric with
%    fluorescent singals
%    FitDataCompare(ShapeMetric, FluorescentMetric)
%      FluorescentMetric is a string containing the name of the fluorescent
%      Metric, eg 'Ff', 'Ffr', 'Ffr.^2', etc.
%   FitDataCompare(ShapeMetric, FluorescentMetric,ShapeBins)

% turn off some warnings
prevWarnState1 = warning('off','MATLAB:load:variableNotFound');
prevWarnState2 = warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');

% here are some default options
if nargin>=4
    options=varargin{4};
else
    options = struct();
end

options = applyDefaults(options);

switch nargin
    case 0
        geometricPropertyName='gc';
        fluorescentMeasurementName='Ff';
        folderPath=uigetdir('/Users/jeffnguyen/Documents/Data/ImagesData/');
        
    case 1
        geometricPropertyName=varargin{1};
        fluorescentMeasurementName='Ff';
        folderPath=uigetdir('/Users/jeffnguyen/Documents/Data/ImagesData/');
        
    case 2
        geometricPropertyName=varargin{1};
        fluorescentMeasurementName=varargin{2};
        folderPath=uigetdir('/Users/jeffnguyen/Documents/Data/ImagesData/');
        
    case {3,4}
        geometricPropertyName=varargin{1};
        fluorescentMeasurementName=varargin{2};
        folderPath=varargin{3};
end

% store the name of the geometric property being used
if isfield(options,'geometricPropertyName')
    if not(isequal(options.geometricPropertyName,geometricPropertyName))
       error('FitDataCompareTri:geometricPropertyMismatch',['The input geometric property in options, ''',...
           options.geometricPropertyName,''' , does not match the explicit input of ''', geometricPropertyName,'''. Aborting']);
    end
end
options.geometricPropertyName = geometricPropertyName;

% store the name of the fluorescent property being used
if isfield(options,'fluorescentMeasurementName')
    if not(isequal(options.fluorescentMeasurementName,fluorescentMeasurementName))
       error('FitDataCompareTri:fluorescentPropertyMismatch',['The input fluorescent property in options, ''',...
           options.fluorescentMeasurementName,''' , does not match the explicit input of ''', fluorescentMeasurementName,'''. Aborting']);
    end
end
options.fluorescentMeasurementName = fluorescentMeasurementName;


undersample=options.undersample;
midflag=options.midflag;
pole_flag=options.poleflag;
plotflag=options.plotflag;
proteinString = options.proteinString;



d=dir(folderPath);
l=[];
NN=[];
N2=[];


for i=1:size(d,1)
    
    % Generate a list of file indices that don't match
    if isempty(regexp(d(i,1).name,'([T][R][I])([ ][\(][0-9]*[\)][.][m][a][t]|[.][m][a][t])','once'))
        l=[l,i];
    end
    
    
end
d(l)=[];
g=[];
gg=[];
aScale=[];
A=[];

% pull out the edges of bins to use for histogram purposes
edgesHist = options.edgesHist;
if isempty(edgesHist)
    % default bins edges
    switch geometricPropertyName
        case 'gc'
%   even if poles are ignored, they should not be removed from the data
%   until after the normalization step
%             if pole_flag==0
%                 edgesHist=-2:1*10^-1:4; %gaussian curvature
%                 
%             else
                edgesHist=-50:2*10^-1:50; %gaussian curvature
                %  edgesHist=-10:2*10^-1:15; %gaussian curvature
                 
%             end
        case {'kminus','k-'}
%             if pole_flag==0
%                 edgesHist=-1:.5*10^-1:1; %k-
%             else
                edgesHist=-1:.5*10^-1:6; %k-
%             end
        case {'kplus','k+'}
            
            edgesHist=-1:1*10^-1:10; %k+
        case 'gm'
%             if pole_flag==0
%                 edgesHist=-0:5*10^-2:2; %mean curvature
%             else
                edgesHist=-0:5*10^-2:10; %mean curvature
%             end
        case 'CLK'
            edgesHist=-0.0001:0.00003:.0025;
        case 'CLT'
            edgesHist=-0.15:.0002:.15;
            edgesHist=-0.01:.0002:.01;
            
        case 'Other'
            edgesHist=0:.05:.999;
        case 'theta'
            edgesHist = -180:2:180;
    end
end

%end
% remove data associated with the > final bin edge
counts=zeros(length(edgesHist')-1,1);

for i=1:size(d,1)
    
    %load data
    try
    cellData=load([folderPath filesep d(i).name],'F','Fcoord','Fcurvature','Ff','Ffr','background');
    catch ME
        cellData=load([folderPath filesep d(i).name],'F','Fcoord','Fcurvature','Ff','Ffr');
    end
    
    try
        F=cellData.F;
        Fcoord=cellData.Fcoord;
        Fcurvature=cellData.Fcurvature;
        Ffr=cellData.Ffr;
        Ff=cellData.Ff;
        
        
    catch ME
        warning('FitDataCompare:fieldsMissing',['Fields are missing from ',d(i).name, '. Skipping and proceeding with next file.']);
        continue
    end
    
    if isfield(cellData,'background')
        options.backgroundValue = cellData.background;
    elseif options.backgroundSubtract == -990;

        warning('FitDataCompare:backgroundMissing',['Background value is missing from ',d(i).name, '. Skipping and proceeding with next file.']);
        continue    
    end
    
    xf=F.vertices(:,1);
    yf=F.vertices(:,2);
    zf=F.vertices(:,3);
    CL=Fcoord.CL;
    gc=Fcurvature.kg; %guassian curv
    gm=Fcurvature.km; %mean curv
    kplus=Fcurvature.k1;
    kminus=Fcurvature.k2;
    
    theta = Fcoord.theta;
    %         Ff=cellF.Ff;
    %         Ffr=cellF.Ffr;
    
    if ~strcmp(fluorescentMeasurementName,'none')
        switch fluorescentMeasurementName
            case 'gc'
                kk = gc;
            case 'amountNotConcentration'
                [~,a,~]=triArea(F);
                
                kk = Ff;
                kk = kk-nanmin(kk(:));
                kk=kk/nanmean(kk(:));
                kk = kk.*a;
            otherwise
                % BPB NOTES: there should be a scaling by area during
                % normalization, but the triangular faces all have similar
                % area
                eval(sprintf('%s=%s;','kk',fluorescentMeasurementName));
                
                switch options.backgroundSubtract
                    % background subtraction options are
                    % -999, use min to max of data
                    % scalar between -100 and +65K, subtract that value from everything
                    % -990, use options.backgroundValue
                    % -998, use min to max of data, normalize by surface area
                    % -995, use nanmin and nanmean to normalize the data
                    case -999
                        kk=kk-nanmin(kk(:));
                        kk=kk/nanmean(kk(:));
                        [~,a,~]=triArea(F);

                    case -998
                        kk = kk-nanmin(kk(:));
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
                    case -995
                        kk=kk-nanmin(kk(:));
                        kk=kk/nanmean(kk(:));
                    case -990
                        kk = kk-options.backgroundValue;
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
                    otherwise
                        kk = kk-options.backgroundSubtract;
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
                end
        end
    else
        
        kk=ones(size(xf));
    end
    
    
    
    if strfind(geometricPropertyName,'CL')
        switch options.backgroundSubtract
                    % background subtraction options are
                    % -999, use min to max of data
                    % scalar between -100 and +65K, subtract that value from everything
                    % -990, use options.backgroundValue
                    % -998, use min to max of data, normalize by surface area
                    % -995, use nanmin and nanmean to normalize the data
                    case -999
                        kk=kk-nanmin(kk(:));
                        kk=kk/nanmean(kk(:));
                    case -998
                        kk = kk-nanmin(kk(:));
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
                    case -995
                        kk=kk-nanmin(kk(:));
                        kk=kk/nanmean(kk(:));
                    case -990
                        kk = kk-options.backgroundValue;
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
                    otherwise
                        kk = kk-options.backgroundSubtract;
                        [~,a,~]=triArea(F);
                        kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
        end
                
        
        
        %% if looking at a centerline property, we'll need to find the
        % centerline index for each surface pixel
        
        % first find the contour length along the CL
        distAlongCL = Fcoord.CL;
        % take the difference along the curve
        distAlongCL = diff(distAlongCL);
        distAlongCL = sqrt(nansum(distAlongCL.^2,2));
        % add up all the distances
        distAlongCL = cumsum(distAlongCL);
        % start at 0
        distAlongCL = [0;distAlongCL];
        
        % for each surface pixel, find the CL point that matches
        [Lia,Lib] = ismember(Fcoord.L,distAlongCL);
        if any(~Lia)
            warning('FitDataCompare:Tri:clIdxMissing','Some surface points do not have a matching centerline index.');
        end
    end
    
    %%
    switch geometricPropertyName
        case 'gc'
            
            k=triSmooth(F.faces,gc)*10^6;% gaussian curvature in 1/um2
        case {'k-','kminus'};
            
            k=triSmooth(F.faces,kminus)*10^3; %k-
        case {'k+','kplus'};
            
            k=triSmooth(F.faces,kplus)*10^3;
        case 'gm'
            
            k=triSmooth(F.faces,gm)*10^3;
            
        case 'CLT'
            clTorsion = Twist(Fcoord.CL,5);
            % evaluate the curvature of the centerline at the appropriate
            % surface indices
            k = clTorsion(Lib);
            %             k=K;
            k=triSmooth(F.faces,k);
            k=k';
        case 'CLK'
            clCurvature = Curvature(Fcoord.CL,5);
            % evaluate the curvature of the centerline at the appropriate
            % surface indices
            k = clCurvature(Lib);
            k=triSmooth(F.faces,k);
            k=k';
            
            % remove pole data
            if pole_flag~=1
                k(Fcurvature.kg>pole_flag) = nan;
            end
        case 'Other'
            if pole_flag==0
                k=trimDataEnds((size(I,2)-size(x0,2))/2,I);
                
            else
                k=I;
            end
            % normalization
            warning('FitDataCompareTri:dataNormalization','Range of data scaled from 0 to 1');
            k=normalizeRange(k);
            
        case 'theta'
            k = triSmooth(F.faces,theta)*180/(pi);
    end
    %%
    % if one is not looking at a centerline property, the string search of
    % geometric property will not include CL
    if isempty(strfind(geometricPropertyName,'CL'))
        
%         % find the surface area of each vertex
%         [~,a,~]=triArea(F);
%         
%         % BPB notes: here's where the normalization takes place on a
%         % per-file (one cell) basis
%         switch fluorescentMeasurementName
%             case {'gc','amountNotConcentration'}
%                 
%             otherwise
%                 kk=normalizeRange(kk);
%                 kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
%         end
        
        
        
        isToRemoveNanOrRange=[];
        % if one only wants to use the center third of the z coordinates,
        % kill off the small and the large Z coordinates
        if midflag
            Z = normalizeRange(zf);
            isToRemoveNanOrRange=(Z<.35 | Z>.65);
        else
            isToRemoveNanOrRange = zeros(size(zf));
        end
        isToRemoveNanOrRange=(isToRemoveNanOrRange>0);
        
        
        if options.poleflag==0
            % if one is supposed to ignore the poles, determine where they are
            prefs.displayFigure = 0;
            prefs.radiusMultiplier = 1.5;
            
            [poleIdx,isPolarRegionVertex,isPolarRegionFace] = poleFinder(F,Fcoord,prefs);
            
            isToRemoveNanOrRange = or(isToRemoveNanOrRange,isPolarRegionVertex);
        elseif options.poleflag == -1
            % remove the poles before the normalization step. Or, remove
            % the poles and then normalize again
            
            prefs.displayFigure = 0;
            prefs.radiusMultiplier = 1.5;
            
            [poleIdx,isPolarRegionVertex,isPolarRegionFace] = poleFinder(F,Fcoord,prefs);
            
            k(isPolarRegionVertex) = [];
            a(isPolarRegionVertex) = [];
            kk(isPolarRegionVertex) = [];
            
            isToRemoveNanOrRange(isPolarRegionVertex) = [];
            % occassionaly there are no non-polar regions
            if isempty(kk)
                warning('FitDataCompareTri:allNan','Following pole removal, there are no remaining surface elements.');
               continue;
            else
                
                % possible place to change
            switch options.backgroundSubtract
                case -999
                kk=normalizeRange(kk);
                otherwise
            end
            
            kk=kk/nansum(kk(:).*a(:)/nansum(a(:)));
            end
        end
    
        % remove these indices from all the data
        k(isToRemoveNanOrRange)=[];
        a(isToRemoveNanOrRange)=[];
        kk(isToRemoveNanOrRange)=[];
        
        
        
        % reshape the data into the proper format
        nSize=numel(k);
        kk=reshape(kk,1,nSize);
        a=(reshape(a,1,nSize));
        k=reshape(k,1,nSize);
        
        % BPB notes, if one is using undersampling, should this be a random
        % subset of vertices instead of every 4? And why should the spacing
        % be exactly 4?
        if undersample
            k=k(1:4:end);
            kk=kk(1:4:end);
            a=a(1:4:end);
        end
    else
        nSize=numel(k);
        kk=reshape(kk,1,nSize);
        a=ones(1,nSize);
        k=reshape(k,1,nSize);
        k=k(1:4:end);
        kk=kk(1:4:end);
        a=a(1:4:end);
        
    end
    
    %%
    
    A=[A,a];
    g=[g,k];
    gg=[gg,kk];
    
    
    [bins,~,~,locs]=w_histcn(k',a',edgesHist);
    % remove data associated with the > final bin edge
    if length(bins)==length(edgesHist)
        bins = bins(1:end-1);
        warning('FitDataCompare:dataExceedsBins',...
            ['Max value for data is ', num2str(max(k(:))),'. ',...
            'Max edge for histogram is ', num2str(edgesHist(end)),'.']);
        
    elseif length(bins)==(length(edgesHist)-1)
    else
        warning('FitDataCompare:wrongSizeBins','Unexpected size of output from w_histcn');
        break
    end
    
    bins=bins*nansum(a);
    counts=counts+bins;
    X=edgesHist;
    
    for ii=2:length(edgesHist)
        Ksample=(k<X(ii) & k>X(ii-1));
        N2(ii-1,i)=nansum(kk(Ksample).*a(Ksample))/nansum(a(Ksample));
    end
    
%     if or(all(isnan( N2(:,i))),all(isnan(kk(Ksample).*a(Ksample))))
%         warning('fitDataCompareTri:AllNansAvg','Average fluorescence values are NaN');
%         disp([folderPath filesep d(i).name]);
%     end
        
    if any(isnan( Ff));
        warning('fitDataCompareTri:SomeNansFf','Some fluorescence values are NaN');
        disp([folderPath filesep d(i).name]);
    end
        
 
    
    
    
    
end
counts=counts/nansum(counts);
switch fluorescentMeasurementName
    case 'gc'
        signalBins = -10:6*10^-1:25;
    case ''
        
    otherwise
        signalBins = 0:.1:2;
end

[count edge mid loc]=w_histcn([g',gg'],A,edgesHist,signalBins);
ctcontrol=w_histcn(gg',A,signalBins);

N3 = nanmean(N2,2);
N3se = nanstd(N2,1,2)./sqrt(nansum(~isnan(N2),2));
X=edgesHist;
N4 = N2;
if ~strcmp(fluorescentMeasurementName,'none')
    
    clear N2
    for ii=1:50
        Krand=randperm(numel(g));
        for i=2:length(X)
            N2(i-1)=nansum(gg(Krand(g<X(i) & g>X(i-1))).*A(Krand(g<X(i) & g>X(i-1))))/nansum(A(Krand(g<X(i) & g>X(i-1))));
        end
        NN(ii,:)=N2;
    end
    for i=2:length(X)
        N2(i-1)=nansum(gg(g<X(i) & g>X(i-1)).*A(g<X(i)&g>X(i-1)))/nansum((g<X(i)&g>X(i-1)));
    end
    N3(N3se==0 | N3se>.5)=NaN;
end

if nargout==0 || plotflag
    
    subplot(2,1,1)
    plot(X(2:end),(mean(NN)+1.96*std(NN)),'r');
    hold on
    plot(X(2:end),(mean(NN)-1.96*std(NN)),'r');
    
    
    xlim([min(X),max(X)]);
    errorbar(X(2:end),(N3),N3se,'.g')
    hold on
    
    ylim([0,2]);
    xlim([min(X),max(X)]);
    % remove data associated with the > final bin edge
    X(end) = [];
    switch geometricPropertyName
        case 'gc'
            xlabel('Gaussian Curvature (1/\mum^2)','FontSize',20);
            ylabel(['Enrichment of ' ,proteinString, ' signal'],'FontSize',14);
            
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Gaussian Curvature (1/\mum^2)','FontSize',20);
        case 'kminus'
            xlabel('kminus (1/nm)');
            ylabel(['Enrichment of ' ,proteinString, ' signal']);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Gaussian Curvature (1/nm^2)');
        case 'kplus'
            xlabel('kplus (1/nm)');
            ylabel(['Enrichment of ' ,proteinString, ' signal']);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('kplus (1/nm)');
        case 'gm'
            xlabel('Mean Curvature (1/\mum)','FontSize',20);
            ylabel(['Enrichment of ' ,proteinString, ' signal'],'FontSize',14);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Mean Curvature (1/\mum)','FontSize',20);
        case 'CLK'
            xlabel('Centerline Curvature(1/nm)','FontSize',14);
            ylabel(['Enrichment of ' ,proteinString, ' signal'],'FontSize',14);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Centerline Curvature(1/nm)','FontSize',14);
        case 'CLT'
            xlabel('Centerline Torsion(1/nm)','FontSize',14);
            ylabel(['Enrichment of ' ,proteinString, ' signal'],'FontSize',14);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Centerline Torsion(1/nm)','FontSize',14);
        case 'Other'
            ylabel(['Enrichment of ' ,proteinString, ' signal']);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Gaussian Curvature (1/nm^2)');
        case theta
            ylabel(['Enrichment of ' ,proteinString, ' signal']);
            subplot(2,1,2);
            bar(X,counts);
            xlim([min(X),max(X)]);
            xlabel('Angular position (degrees)');
    end
    
    ylabel('Frequency','FontSize',14);
    
else
end

% remove nans and data outside the range of the bins
isToRemoveNanOrRange=(isnan(g)|A==0|g>max(X)|g<min(X));
A(isToRemoveNanOrRange) = [];
g(isToRemoveNanOrRange) = [];
gg(isToRemoveNanOrRange) = [];

%% save data

if length(X) > length(N3);
    N3 = [N3;nan(1,length(X)-length(N3))];
end
if length(X) > length(counts);
    counts = [counts;nan(1,length(X)-length(counts))];
    
end

        
Fittingdata.N3 = N3;
Fittingdata.counts = counts;
Fittingdata.X = X';
Fittingdata.N3se = N3se;
[parentFolder,strain] = fileparts(folderPath);
[parentFolder,dateTaken] = fileparts(parentFolder);

nColumns = any(not(isnan(N4)));
nCells = sum(nColumns);

Fittingdata.label = [dateTaken,'-',strain,'n=',num2str(nCells)];

modeName = enrichmentModeString(options);
output = [parentFolder,filesep,dateTaken,'-',strain,'_Mode-',modeName,'_curve.mat'];
save(output,'Fittingdata')


warning(prevWarnState1);
warning(prevWarnState2);

    

function options = applyDefaults(options)

if not(isfield(options,'backgroundSubtract'))
    % background subtraction options are 
    % -999, use min to max of data
    % scalar between -100 and +65K, subtract that value from everything
    % -990, use options.backgroundValue
    % -998, use min to max of data, normalize by surface area
    % -995, use nanmin and nanmean to normalize the data
    options.backgroundSubtract = -999;
end


if not(isfield(options,'undersample'))
    % undersample 0 uses everything, undersample 1 uses 1/4 of the points
    % don't use every surface point as there is still some level of correlation in surface intensity values
    options.undersample=0;
end

if not(isfield(options,'midflag'))
    % mid flag 1 uses data only from the middle Z sections
    options.midflag=0;
end


if not(isfield(options,'poleflag'))
    % pole flag 1 includes the poles in the data, 0 discards them after the
    % normalization step, -1 discards them before normalization
    options.poleflag=1;
end

if not(isfield(options,'plotflag'))
    % plot flag 1 makes simple plots of enrichment of some signal as a function
    % of some geometric parameter and the distribution of that parameter
    options.plotflag=1;
end

if not(isfield(options,'proteinString'))
    % default string for the name of what the fluorescent protein/stain is
    options.proteinString='proteinNameHere';
end

if not(isfield(options,'edgesHist'))
    % default to empty bins
    options.edgesHist = [];
end
