function [iSeries,nSeries,iPlane,nPlanes,iColor,nColors,iTime,nTimes,stackString] = parseND2label(nd2label);
%%
nd2label = textscan(nd2label,'%s','delimiter',';');
nd2label = nd2label{1};
%% series : index 2
seriesString = nd2label{2};
seriesString = textscan(seriesString,'%s');
seriesString = seriesString{1};
seriesString = textscan(seriesString{2},'%f','delimiter','/');
iSeries = seriesString{1}(1);
nSeries = seriesString{1}(2);
%% plane : index 3
planeString = nd2label{3};
planeString = textscan(planeString,'%s');
planeString = planeString{1};
planeString = textscan(planeString{2},'%f','delimiter','/');
iPlane = planeString{1}(1);
nPlanes = planeString{1}(2);
%% color : index 4
colorString = nd2label{4};
colorString = textscan(colorString,'%s','delimiter','=');
colorString = colorString{1};
colorString = textscan(colorString{2},'%f','delimiter','/');
iColor = colorString{1}(1);
nColors = colorString{1}(2);
%% time : index 5
timeString = nd2label{5};
timeString = textscan(timeString,'%s','delimiter','=');
timeString = timeString{1};
timeString = textscan(timeString{2},'%f','delimiter','/');
iTime = timeString{1}(1);
nTimes = timeString{1}(2);
%% composite string
fileLabel = nd2label{1};
[~,filename] = fileparts(fileLabel);
stackString = [filename,'_series',num2str(iSeries),'_C=',num2str(iColor)];
