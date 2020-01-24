function displayPlots_signalSingleCells(summaryCellArray,nSplits,maxSignalTotal)
%%
% how to estimate the appropriate bin width for histograms?
nFolders = size(summaryCellArray,1);
for iiFolder  = 1:nFolders
   nFiles(iiFolder) = size(summaryCellArray{iiFolder,3},1); 
end

%  K = 1 + log2(nFiles);
% K = 2*nFiles.^(1/3);
% nBins = round(mean(K));

nBins = 100;

% data can sometimes be on ridiculously different scales, leading to silly
% looking plots



% display a kernel smoothing histogram
figure();
for iiChannel = 1:nSplits
    subplot(1,nSplits,iiChannel);
hold on;
   binEdges = linspace(1,log2(1.25*maxSignalTotal(iiChannel)),nBins);
   for iiFolder = 1:nFolders
   densityTemp = ksdensity(log2(summaryCellArray{iiFolder,3}(:,iiChannel)),binEdges);
   % normalize to probability from probability density
   densityTemp = densityTemp./sum(densityTemp);
   plot(2.^binEdges,densityTemp,'displayname',summaryCellArray{iiFolder,4});
   end
end
        set(gca,'xscale','log');
legend();

figure();
for iiChannel = 1:nSplits
    axesHand = subplot(2,nSplits,iiChannel);
    
    for iiFolder = 1:nFolders
        tempCellArrayData{iiFolder} = summaryCellArray{iiFolder,3}(:,iiChannel);
        tempCellLabels{iiFolder} = summaryCellArray{iiFolder,4};
    end
    
   violinHandle = violin(tempCellArrayData);
   if verLessThan('MATLAB','R2018b');
       set(axesHand,'XTick',[1:nFolders]);
       set(axesHand,'XTickLabels',tempCellLabels);
       set(axesHand,'XTickLabelRotation',45);
       set(axesHand,'TickLabelInterpreter','none');

   else
       xticks(axesHand,[1:nFolders]);
       xticklabels(axesHand,tempCellLabels);
       xtickangle(axesHand,45);
       set(axesHand,'TickLabelInterpreter','none');
   end
   
  childHandles = get(axesHand,'children');
  minVal = inf;
   for ii = 1:size(childHandles,1)
       switch class(childHandles(ii))
           case 'matlab.graphics.chart.primitive.Line'
               tempYdata = get(childHandles(ii),'YData');
               tempYdata(tempYdata<=0) = min(tempYdata>0);
               minVal = min(minVal,min(tempYdata>0));
               set(childHandles(ii),'YData',tempYdata);
           case 'matlab.graphics.primitive.Patch'
               tempYdata = get(childHandles(ii),'Vertices');
               replaceIndex = tempYdata(:,2)<=0;
               tempYdata(replaceIndex,2) = min(tempYdata(not(replaceIndex),2));
               set(childHandles(ii),'Vertices',tempYdata);
               minVal = min(minVal, min(tempYdata(not(replaceIndex),2)));
           otherwise
               keyboard
       end
%        tempYdata = get(violinHandle(ii),'Vertices');
%        replaceIndices = tempYdata(:,2)<=0;
%        tempYdata(replaceIndices,2) = min(tempYdata(not(replaceIndices),2));
%        set(violinHandle(ii),'Vertices',tempYdata);
   end
   
   

   set(axesHand,'Yscale','log');
%    ylimA = ylim(axesHand);
%    ylim(axesHand,[minVal,ylimA(2)]);
   
end