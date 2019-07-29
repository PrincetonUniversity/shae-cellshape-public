function  plotter2(hObject,eventdata)
% get the workspace ready
clear Ff Ffr

% switch the loading variable not found warning to an error
warnstatus = warning();
warning('error','MATLAB:load:variableNotFound');
warning('off','MATLAB:Axes:UpVector');
%% find which gui this came from and the handles structure of that gui
handles=guidata(get(hObject,'Parent'));
currentImageIdx=getappdata(handles.figure1,'currentImageIdx');
imageList=getappdata(handles.figure1,'imageList');
fileNameToDisplay = imageList(currentImageIdx).name;
folder=getappdata(handles.figure1,'folder');

setappdata(handles.figure1,'currentImage',imageList(currentImageIdx).name);
fullFileName = fullfile(folder,fileNameToDisplay);
setappdata(hObject,'currentImage',fullFileName);

% we should only load the pieces of data that we need. This should improve
% the speed of the load

% first, check and see which data type we're dealing with
tail = get(handles.fileEnder,'String');
displayStrings = get(handles.popupmenu1,'String');
displayType = displayStrings{get(handles.popupmenu1,'Value')};
defaultTypeIndx = find(strcmp(displayStrings,'Wireframe only'));

switch tail
    case 'TRI'
        % if loading in the triangular meshes, check and see if we need the
        % fluorescent signal or not
        switch displayType
            case {'Ff'};
                try
                    load(fullFileName,'Fcoord','F','Ff');
                catch ME
                    msgbox(['Unable to use ''',displayType,'''. Switching to ''wireframe only''']);
                    load(fullFileName,'Fcoord','F');
                    set(handles.popupmenu1,'Value',defaultTypeIndx)
                    displayType = displayStrings{defaultTypeIndx};
                end
            case {'Ffr'};
                try
                    load(fullFileName,'Fcoord','F','Ffr');
                catch ME
                    msgbox(['Unable to use ''',displayType,'''. Switching to ''wireframe only''']);
                    load(fullFileName,'Fcoord','F');
                    set(handles.popupmenu1,'Value',defaultTypeIndx)
                    displayType = displayStrings{defaultTypeIndx};
                end
            case {'Gaussian curvature'}
                load(fullFileName,'Fcoord','F','Fcurvature');
            otherwise
                load(fullFileName,'Fcoord','F');
                
        end
                
                
    otherwise
        % go ahead and load everything if it's not clear what to do
        
        load(fullFileName)
        
        if exist('cellcoord','var');
            xf=cellcoord.x;
            yf=cellcoord.y;
            zf=cellcoord.z;
            CL=cellcoord.CL;
            gc=cellcurve.kg;
            gm=cellcurve.km;
            kplus=cellcurve.k1;
            kminus=cellcurve.k2;
            r=cellShape.Radius;
            clear Ff Ffr
            try
                Ff=cellF.Ff2;
                Ffr=cellF.Ffr2;
                cellF.mFraction;
                if size(Ff,1)==size(xf,1)-1
                    Ff=[Ff(end,:);Ff];
                end
                if size(Ffr,1)==size(xf,1)-1
                    Ffr=[Ffr(end,:);Ffr];
                end
            catch ME
            end
            
            % if the interpolated intensities don't exist, use the z position
            if not(exist('Ff','var'))
                Ff = zf;
            end
            
            if not(exist('Ffr','var'))
                Ffr = zf;
            end
            
            
        else
            polyFit=[];
        end
        
        %% this is the end of the otherwise
end

% once again check and see which mode is supposed to be displayed. If there
% had been a loading error, the mode may have changed.
displayStrings = get(handles.popupmenu1,'String');
displayType = displayStrings{get(handles.popupmenu1,'Value')};
defaultTypeIndx = find(strcmp(displayStrings,'Wireframe only'));

%%
% determine which display mode we should be using
if exist('Fcoord','var');
    displayMode = 'patch';
else
    displayMode = 'surface';
end

%% mean center the CL and the surface
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
        % why might there be no centerline?
        if isfield(Fcoord,'CL')
            Fcoord.CL = bsxfun(@minus,Fcoord.CL,mean(F.vertices));
        else
            warning('FitViewer:plotter2:no_CL',['No centerline in file: ',fullFileName]);
        end
        F.vertices = bsxfun(@minus,F.vertices,mean(F.vertices));
        
end






%% display the surface
switch displayMode
    case 'surface'
        %  switch col
        %      case 'Z'
        
        fig = surf(hObject,xf,yf,zf,Ffr,'EdgeColor','none');
        
        
        
        
        %      case 'Gaussian Curvature'
        %          fig=surf(handles.axes1,xf,yf,zf, gc,'EdgeColor','none');
        %      case 'Mean Curvature'
        %          fig=surf(handles.axes1,xf,yf,zf,gm,'EdgeColor','none');
        %      case 'Fluorescent Signal'
        %          %     sFf=wiener2(repmat(imtophat(Ff,se),3,1));
        %          %     sFf=sFf(size(Ff,1)+1:2*size(Ff,1),:);
        %          fig=surf(handles.axes1,xf,yf,zf,Ffr,'EdgeColor','none');
        %          %surf(handles.axes1,sFf,'EdgeColor','none');
        %      case 'K plus'
        %          fig=surf(handles.axes1,xf,yf,zf,kplus,'EdgeColor','none');
        %      case 'K minus'
        %          fig=surf(handles.axes1,xf,yf,zf,kminus,'EdgeColor','none');
        %      case 'F unwrap'
        %          fig=imagesc(Ff, 'Parent',handles.axes1);
        %      case 'Radius'
        %   fig=surf(handles.axes1,xf,yf,zf,r,'EdgeColor','none');
        %
        %      case 'CL Curvature'
        %       %   fig=surf(handles.axes1,xf,yf,zf,repmat(smooth(K,4)',[size(xf,1),1]),'EdgeColor','none');
        %
        %          fig=tubeplot(CL',20,20,1,K,handles.axes1);
        %      case 'CL Twist'
        %         % fig=surf(handles.axes1,xf,yf,zf,double(repmat(T,[size(xf,1),1])),'EdgeColor','none');
        %
        %              fig=tubeplot(CL',20,20,1,T',handles.axes1);
        %
        %      case 'CL Mean F'
        %              fig=tubeplot(CL',20,20,1,mean(Ffr)',handles.axes1);
        %
        %      otherwise
        %
        %     fig=surf(handles.axes1,xf,yf,zf,'EdgeColor','none');
        % end
        %
        % if ~strcmp(col,'F unwrap')
        hold(hObject,'on')
        
        if ~isempty(polyFit);
            for pi=1:length(polyFit);
                plot3(hObject,polyFit(pi).xs-cellcm(1),polyFit(pi).ys-cellcm(2),...
                    polyFit(pi).zs-cellcm(3),'black','LineWidth',3);
            end
            
        end
        
        % BPB notes, this is defined for rectangular grids of size 40 or 41, should
        % depend on the size of the grid
        loopSize = size(xf,1)-1;
        piOverTwo= round(loopSize/4);
        plot3(hObject,xf(1,:),yf(1,:),zf(1,:))
        plot3(hObject,xf(1+piOverTwo,:),yf(1+piOverTwo,:),zf(1+piOverTwo,:));
        plot3(hObject,xf(1+2*piOverTwo,:),yf(1+2*piOverTwo,:),zf(1+2*piOverTwo,:));
        plot3(hObject,xf(1+3*piOverTwo,:),yf(1+3*piOverTwo,:),zf(1+3*piOverTwo,:));
        axis(hObject,'equal')
        
        hold(hObject,'off')
        
        
        
        % end
        %
        % if get(handles.Colourbar,'Value')
        %     colorbar('peer',handles.axes1);
        % end
        % minlim=min([xlim(handles.axes1),ylim(handles.axes1)]);
        % maxlim=max([xlim(handles.axes1),ylim(handles.axes1)]);
        %
        % xlim(handles.axes1,[minlim maxlim]);
        % ylim(handles.axes1,[minlim maxlim]);
        %
        %
        % % camlight(handles.axes1,'right');
        % % camlight(handles.axes1,'left');
        % % set(fig,'FaceLighting','phong','FaceColor','interp',...
        % %       'AmbientStrength',0.5)
        % %  set(camlight('right'),'Parent',handles.axes1);
        % % set(camlight('left'),'Parent',handles.axes1);
        %
        %
        % if getappdata(handles.figure1,'holdaxes')
        % view(handles.axes1,[az,el])
        %     caxis(handles.axes1,c)
        %   camtarget(handles.axes1,Tcam)
        %  camva(handles.axes1,V)
        % axis(handles.axes1,v)
        %
        % end
        % %set(handles.axes1,'EdgeColor','b')
        % %zoom(handles.axes1,z);
        %
        % setappdata(handles.figure1,'CurrentFrame',round(k));
        % set(handles.jumpCurrentFrame,'string',num2str(round(k)));
        %
        %
        % getappdata(handles.figure1,'CurrentFrame');
        
    case 'patch'
        % %         fig=surf(hObject,xf,yf,zf,'EdgeColor','none');
        % clear the axis
        cla(hObject);
        hold(hObject,'on');
        switch displayType
            case 'Z'
                h = patch(F,'facevertexcdata',F.vertices(:,3),'facealpha',1,'facecolor','interp','edgecolor','none','Parent',hObject);
            case 'Ffr'
                h = patch(F,'facevertexcdata',Ffr,'facealpha',1,'facecolor','interp','edgecolor','none','Parent',hObject);
            case 'Ff'
                h = patch(F,'facevertexcdata',Ff,'facealpha',1,'facecolor','interp','edgecolor','none','Parent',hObject);
            case 'Gaussian curvature'
                h = patch(F,'facevertexcdata',Fcurvature.kg,'facealpha',1,'facecolor','interp','edgecolor','none','Parent',hObject);
            otherwise
                h = patch(F,'facecolor','w','facealpha',.8,'Parent',hObject);
        end
        %        %   h=patch(F2,'facevertexCdata',vertAreaI,'facecolor','interp','facealpha',.5,'edgecolor','none');  axis equal;  view(3);
        axis(hObject,'equal');
        view(hObject,3);
        % no need to flush the drawing queue, it should update fast enough
        
        if get(handles.checkbox2,'Value')==get(handles.checkbox2,'Max')
        drawnow;
        end
        hold(hObject,'off');
        
end


% increment the file index to move on the next file
setappdata(handles.figure1,'currentImageIdx',currentImageIdx+1);
warning(warnstatus);