function  plotter(hObject,eventdata)

clear Ff Ffr
handles=guidata(get(hObject,'Parent'));
k=round(get(hObject,'Value'));
p=getappdata(handles.figure1,'folderlist');
% get a directory listing of the files matching this search term
d=dir([p,filesep,'*',getappdata(handles.figure1,'fileEnder'),'*']);

setappdata(handles.figure1,'totalFrames',length(d));
set(handles.edit3,'String',num2str(length(d)));
k=min(k,length(d));
load([p filesep d(k).name])
if exist('cellF','var');
    xf=cellcoord.x;
    yf=cellcoord.y;
    zf=cellcoord.z;
    CL=cellcoord.CL;
    gc=cellcurve.kg;
    gm=cellcurve.km;
    kplus=cellcurve.k1;
    kminus=cellcurve.k2;
    r=cellShape.Radius;
    if ~isempty(cellF)
        Ff=cellF.Ff2;
        Ffr=cellF.Ffr2;
        if size(Ff,1)==size(xf,1)-1
            Ff=[Ff(end,:);Ff];
        end
        if size(Ffr,1)==size(xf,1)-1
            Ffr=[Ffr(end,:);Ffr];
        end
    end
    
    
else
    polyFit=[];
end



set(handles.imageFile,'String',d(k).name);

if exist('Fcoord','var')
    displayMode = 'patch';
    
else
    displayMode = 'surface';
end

col=getappdata(handles.figure1,'cmap');
v=axis(handles.axes1);
c=caxis(handles.axes1);
V=camva(handles.axes1);
Tcam=camtarget(handles.axes1);
[az,el] = view(handles.axes1);
se = strel('disk',5);

% trim data (not yet updated for patches)
trim=str2num(get(handles.trim,'String'));
if not(isnan(trim))
    switch displayMode
        case 'surface'
            trim=min(round(size(xf,2)/2)-5,trim);
            
            % check to see if there is a fluorescence channel
            if ~isempty(trim) && exist('Ff','var')
                [xf,yf,zf,Ff,Ffr,gc,gm,kplus,kminus,r]=trimDataEnds(trim,xf,yf,zf,Ff,Ffr,gc,gm,kplus,kminus,r);
            elseif ~isempty(trim)
                [xf,yf,zf,gc,gm,kplus,kminus,r]=trimDataEnds(trim,xf,yf,zf,gc,gm,kplus,kminus,r);
                
            end
            
        case 'patch'
            warning('Fit_Viewer:plotter:trim:patch',...
                'trimming Data for patch objects is not yet supported')
    end
    
end

if get(handles.rowNormalize,'Value')
    Ff=bsxfun(@rdivide,Ff,median(Ff,2));
    Ffr=bsxfun(@rdivide,Ff,median(Ffr,2));
    
end

% mean center
switch displayMode
    case 'surface'
        cellcm=[mean(xf(:)),mean(yf(:)),mean(zf(:))];
        if exist('CL','var')
            CL(:,1)=CL(:,1)-cellcm(1);
            CL(:,2)=CL(:,2)-cellcm(2);
            CL(:,3)=CL(:,3)-cellcm(3);
        end
        xf=xf-cellcm(1);
        yf=yf-cellcm(2);
        zf=zf-cellcm(3);
        
    case 'patch'
        if not(isfield(Fcoord,'CL'));
            warning('Fit_Viewer:plotter:noCL',...
                ['No CL for file: ' p filesep d(k).name]);
        else
            Fcoord.CL = bsxfun(@minus,Fcoord.CL,mean(F.vertices,1));
        end
        F.vertices = bsxfun(@minus, F.vertices,mean(F.vertices,1));
        
end
fig=getappdata(0,'fig');
if ~isempty(fig) && all(ishandle(fig)), delete(fig); end

% remove the patch if it exists
switch displayMode
    case 'surface'
    case 'patch'
        cla(handles.axes1);
        %         patchHandle = findobj(handles.axes1,'Type','patch');
        %         if not(isempty(patchHandle))
        %             delete(patchHandle);
        %         end
        
        %         % tube plots?
        %         tubeHandle = findobj(handles.axes1,'Tag','tubePlot');
        %         if not(isempty(tubeHandle))
        %             delete(tubeHandle);
        %         end
        
end


% each type of plotting is slightly different
try
    switch col
        case 'Z'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf,'EdgeColor','none');
                case 'patch'
                    patch(F,'facevertexCdata',F.vertices(:,3),'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
            
        case 'Gaussian Curvature'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf, gc,'EdgeColor','none');
                case 'patch'
                    patch(F,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
            
        case 'Mean Curvature'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf,gm,'EdgeColor','none');
                    
                case 'patch'
                    patch(F,'facevertexCdata',Fcurvature.km,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
            
        case 'Fluorescent Signal'
            switch displayMode
                case 'surface'
                    %     sFf=wiener2(repmat(imtophat(Ff,se),3,1));
                    %     sFf=sFf(size(Ff,1)+1:2*size(Ff,1),:);
                    fig=surf(handles.axes1,xf,yf,zf,Ff,'EdgeColor','none','facecolor','inter','facealpha',.75);
                    xlabel(handles.axes1,'x');
                    ylabel(handles.axes1,'y');
                    zlabel(handles.axes1,'z');
                    %surf(handles.axes1,sFf,'EdgeColor','none');
                case 'patch'
                    patch(F,'facevertexCdata',Ff,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
            if ~isempty(polyFit);
                hold(handles.axes1,'on');
                for pi=1:length(polyFit);
                    plot3(handles.axes1,polyFit(pi).xs-cellcm(1),polyFit(pi).ys-cellcm(2),...
                        polyFit(pi).zs-cellcm(3),'black','LineWidth',3);
                end
                hold(handles.axes1,'off');
                
            end
        case 'Fluorescent Signal Rprojected'
            switch displayMode
                case 'surface'
                    %     sFf=wiener2(repmat(imtophat(Ff,se),3,1));
                    %     sFf=sFf(size(Ff,1)+1:2*size(Ff,1),:);
                    fig=surf(handles.axes1,xf,yf,zf,Ffr,'EdgeColor','none','facecolor','inter','facealpha',.75);
                    xlabel(handles.axes1,'x');
                    ylabel(handles.axes1,'y');
                    zlabel(handles.axes1,'z');
                    %surf(handles.axes1,sFf,'EdgeColor','none');
                case 'patch'
                    fig=patch(F,'facevertexCdata',Ffr,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
            if ~isempty(polyFit);
                hold(handles.axes1,'on');
                for pi=1:length(polyFit);
                    plot3(handles.axes1,polyFit(pi).xs-cellcm(1),polyFit(pi).ys-cellcm(2),...
                        polyFit(pi).zs-cellcm(3),'black','LineWidth',3);
                end
                hold(handles.axes1,'off');
                
            end
            
        case 'K plus'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf,kplus,'EdgeColor','none');
                case 'patch'
                    fig=patch(F,'facevertexCdata',Fcurvature.k2,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
            
        case 'K minus'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf,kminus,'EdgeColor','none');
                    
                case 'patch'
                    fig=patch(F,'facevertexCdata',Fcurvature.k1,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
            
        case 'F unwrap'
            switch displayMode
                case 'surface'
                    end_pts=size(xf,2)-size(x0,2);
                    end_pts=round(end_pts/2);
                    [Ff0]=trimDataEnds(end_pts,Ff);
                    fig=imagesc(Ff0, 'Parent',handles.axes1);
                    set(handles.axes1,'ydir','normal');
                case 'patch'
                    warning('Fit_Viewer:plotter:patch:Unwrap',...
                        'Unwrap images are not yet available for patch objects');
                    
                    %                 patch(F,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
            
        case 'F thresh'
            switch displayMode
                case 'surface'
                    %                 end_pts=size(xf,2)-size(x0,2);
                    %                 end_pts=round(end_pts/2);
                    %                 [x0,y0,z0,Ff0]=...
                    %                     trimDataEnds(end_pts,xf,yf,zf,Ff);
                    if exist('spotFind4','file')
                        
                        FZ=spotFind4(Ffr,1);
                        
                        FZ=max(circshift(FZ,size(FZ,1)/3),FZ);
                        FZ=max(circshift(FZ,size(FZ,1)/3),FZ);
                        FZ=FZ(1:size(Ff,1),:);
                        % FZ=bpass_jn(Ff,.5,[20,20]);
                        %                      FZ=eigenSegmentation(FZ);
                        fig=surf(handles.axes1,xf,yf,zf,(FZ>0).*Ff,'EdgeColor','none');
                        if ~isempty(polyFit);
                            hold(handles.axes1,'on');
                            for pi=1:length(polyFit);
                                plot3(handles.axes1,polyFit(pi).xs-cellcm(1),polyFit(pi).ys-cellcm(2),...
                                    polyFit(pi).zs-cellcm(3),'black','LineWidth',3);
                            end
                            hold(handles.axes1,'off');
                            axis(handles.axes1,'equal','tight')
                            
                        end
                        
                    else
                        fig=imagesc(Ff0, 'Parent',handles.axes1);
                        warning('Fit_Viewer:plotter:threshold','You do not have the Thresholding Function on your path');
                    end
                    
                    set(handles.axes1,'ydir','normal');
                case 'patch'
                    warning('Fit_Viewer:plotter:patch:Unwrap',...
                        'Unwrap images are not yet available for patch objects');
                    
                    %                 patch(F,'facevertexCdata',Fcurvature.kg,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
            
        case 'Radius'
            switch displayMode
                case 'surface'
                    fig=surf(handles.axes1,xf,yf,zf,r,'EdgeColor','none');
                case 'patch'
                    fig=patch(F,'facevertexCdata',Fcoord.r,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
            end
        case 'Theta'
            switch displayMode
                case 'surface'
                    %     sFf=wiener2(repmat(imtophat(Ff,se),3,1));
                    %     sFf=sFf(size(Ff,1)+1:2*size(Ff,1),:);
                    %                 fig=surf(handles.axes1,xf,yf,zf,Ff,'EdgeColor','none');
                    %                 xlabel(handles.axes1,'x');
                    %                 ylabel(handles.axes1,'y');
                    %                 zlabel(handles.axes1,'z');
                    %surf(handles.axes1,sFf,'EdgeColor','none');
                    warning('Fit_Viewer:plotter:surface:theta',...
                        'Displaying the angular coordinate theta is not yet available in surface mode.');
                case 'patch'
                    fig=patch(F,'facevertexCdata',Fcoord.theta,'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
            
        case 'CL Curvature'
            switch displayMode
                case 'surface'
                    %   fig=surf(handles.axes1,xf,yf,zf,repmat(smooth(K,4)',[size(xf,1),1]),'EdgeColor','none');
                    if exist('cellcurve','var')
                        if isfield(cellcurve,'curvature');
                            tubeplot(CL',20,20,1,cellcurve.curvature,handles.axes1);
                            
                        end
                        
                    else
                        tubeplot(CL',20,20,1,K,handles.axes1);
                        
                    end
                    set(get(handles.axes1,'Children'),'Tag','tubePlot');
                    
                case 'patch'
                    K = double(Curvature(Fcoord.CL,4));
                    tubeplot(Fcoord.CL',20,20,1,K,handles.axes1);
                    set(get(handles.axes1,'Children'),'Tag','tubePlot');
            end
            
        case 'CL Twist'
            switch displayMode
                case 'surface'
                    % fig=surf(handles.axes1,xf,yf,zf,double(repmat(T,[size(xf,1),1])),'EdgeColor','none');
                    
                    if exist('cellcurve','var')
                        if isfield(cellcurve,'twist');
                            tubeplot(CL',20,20,1,cellcurve.twist',handles.axes1);
                            
                        end
                        
                    else
                        tubeplot(CL',20,20,1,T,handles.axes1);
                        
                    end
                    set(get(handles.axes1,'Children'),'Tag','tubePlot');
                    
                case 'patch'
                    % the new Twist calculator does not use a smoothing
                    % lengthscale but rather a smoothing tolerance for
                    % deviations between the observed CL positions and the
                    % smoothing spline. By default, this is 1/3 the distance
                    % between subsequent centerline positions
                    % T = double(Twist(Fcoord.CL,4));
                    T = double(Twist(Fcoord.CL));
                    
                    %% calculating Twist is particularly hard at the poles
                    % because there are end effects for calculating derivatives
                    % as well as the contour is forced to go through the pole
                    CL_contourLength = diff(Fcoord.CL);
                    CL_contourLength = sqrt(sum(CL_contourLength.^2,2));
                    CL_contourLength = cumsum(CL_contourLength);
                    CL_contourLength = [0;CL_contourLength];
                    
                    % find the gaussian curvature at the poles
                    [~,poleIDX1] = min(Fcoord.L);
                    [~,poleIDX2] = max(Fcoord.L);
                    maxPoleRadius = sqrt(1./min(Fcurvature.kg([poleIDX1,poleIDX2])));
                    
                    % find how many points 2X R is
                    try
                        
                        firstPostPole = find(CL_contourLength>2*maxPoleRadius,1, 'first');
                        lastPostPole = find((CL_contourLength(end)-CL_contourLength)>2*maxPoleRadius,1, 'last');
                        T(1:(firstPostPole-1)) = nan;
                        T((lastPostPole+1):end) = nan;
                    catch ME
                        T(:) = 0;
                        warning('plotter:Twist:tooMuchPole','No centerline outside of the ''region'' of the poles');
                    end
                    %%
                    tubeplot(Fcoord.CL',20,20,1,T',handles.axes1);
                    set(get(handles.axes1,'Children'),'Tag','tubePlot');
                    
            end
        case 'CL Mean F'
            %%
            switch displayMode
                case 'surface'
                    tubeplot(CL',20,20,1,mean(Ffr)',handles.axes1);
                case 'patch'
                    
                    
                    % find the list of centerline coordinates
                   CL_Llengths = [0;cumsum(sqrt(sum(diff(Fcoord.CL).^2,2)))];
                   % there is some rounding error so every point doesn't map exactly to a centerline coordinate
                   for jjSurfInd = 1:length(Fcoord.L)
                      [~,Fcoord.L_indx(jjSurfInd)] = min(abs(Fcoord.L(jjSurfInd)-CL_Llengths));
                   end
                   
                   % find the surface area of each patch region
                   
                   [~,vertArea,~]=triArea(F);    
                   for jjLIndex  = 1:length(Fcoord.CL)
                      intensityMap(jjLIndex) = ... % average is weighted by surface area at each vertex
                          sum( vertArea(Fcoord.L_indx==jjLIndex).*Ff(Fcoord.L_indx==jjLIndex))...
                          / sum( vertArea(Fcoord.L_indx==jjLIndex));
                   end
                   
                   [~,~]=tubeplot(Fcoord.CL',20,20,1,double(intensityMap'),handles.axes1);

%                     warning('Fit_Viewer:plotter:patch:CLFluor',...
%                         'Centerline fluorescence are not yet available for patch objects');
                    set(get(handles.axes1,'Children'),'Tag','tubePlot');
                    
            end
            
        otherwise
            switch displayMode
                case 'surface'
                    
                    fig=surf(handles.axes1,xf,yf,zf,'EdgeColor','none');
                case 'patch'
                    fig=patch(F,'facevertexCdata',F.vertices(:,3),'facecolor','interp','facealpha',1,'edgecolor','none','Parent',handles.axes1);
                    
            end
    end
    
    if ~strcmp(col,'F unwrap') && ~strcmp(col,'F thresh')
        hold(handles.axes1,'on')
        
        %% plot the "fixed angular coordinate" lines
        switch displayMode
            case 'surface'
                % BPB notes, this is defined for rectangular grids of size 40 or 41, should
                % depend on the size of the grid
                loopSize = size(xf,1)-1;
                piOverTwo= round(loopSize/4);
                if ~strcmp(col,'F unwrap') && ~strcmp(col,'F thresh');
                    plot3(handles.axes1,xf(1,:),yf(1,:),zf(1,:))
                    plot3(handles.axes1,xf(1+piOverTwo,:),yf(1+piOverTwo,:),zf(1+piOverTwo,:));
                    plot3(handles.axes1,xf(1+2*piOverTwo,:),yf(1+2*piOverTwo,:),zf(1+2*piOverTwo,:));
                    plot3(handles.axes1,xf(1+3*piOverTwo,:),yf(1+3*piOverTwo,:),zf(1+3*piOverTwo,:));
                end
            case 'patch'
                % define a window for angles allowed, pi/12 off of pi/2
                pi=3.14159;
                isPointOfInterest = mod(Fcoord.theta,pi/2)<(pi/10);
                scatter3(handles.axes1,F.vertices(isPointOfInterest,1),...
                    F.vertices(isPointOfInterest,2),...
                    F.vertices(isPointOfInterest,3),'b*');
                %         warning('Fit_Viewer:plotter:noFixedAngle',...
                %             'Patch object fixed angular lines have not been enabled');
        end
        axis(handles.axes1,'equal','tight')
        
        hold(handles.axes1,'off')
        setappdata(0,'fig',fig);
    end
    
catch ME2
    if strcmp(col,'Z')
        msgbox(['Unable to use ''Z''. Aborting.']);
        return
    else
        msgbox(['Unable to use ''',col,'''. Switching to ''Z''']);
        % list of fields for display drop down
        listOfStrings = get(handles.popupmenu1,'String');
        newIdx = find(strcmp(listOfStrings,'Z'));
        set(handles.popupmenu1,'Value',newIdx);
        setappdata(handles.figure1,'cmap','Z');
        % rerun the plotting function
        plotter(hObject,eventdata)
        return
    end
end

if get(handles.Colourbar,'Value')
    colorbar('peer',handles.axes1);
end
minlim=min([xlim(handles.axes1),ylim(handles.axes1)]);
maxlim=max([xlim(handles.axes1),ylim(handles.axes1)]);

xlim(handles.axes1,[minlim maxlim]);
ylim(handles.axes1,[minlim maxlim]);


% camlight(handles.axes1,'right');
% camlight(handles.axes1,'left');
% set(fig,'FaceLighting','phong','FaceColor','interp',...
%       'AmbientStrength',0.5)
%  set(camlight('right'),'Parent',handles.axes1);
% set(camlight('left'),'Parent',handles.axes1);


if getappdata(handles.figure1,'holdaxes')
    view(handles.axes1,[az,el])
    caxis(handles.axes1,c)
    camtarget(handles.axes1,Tcam)
    camva(handles.axes1,V)
    axis(handles.axes1,v)
    
end
%set(handles.axes1,'EdgeColor','b')
%zoom(handles.axes1,z);

setappdata(handles.figure1,'CurrentFrame',round(k));
set(handles.jumpCurrentFrame,'string',num2str(round(k)));


getappdata(handles.figure1,'CurrentFrame');