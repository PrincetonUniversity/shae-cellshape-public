function varargout = ScreenFits(varargin)
% SCREENFITS MATLAB code for ScreenFits.fig
%      SCREENFITS, by itself, creates a new SCREENFITS or raises the existing
%      singleton*.
%
%      H = SCREENFITS returns the handle to a new SCREENFITS or the handle to
%      the existing singleton*.
%
%      SCREENFITS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SCREENFITS.M with the given input arguments.
%
%      SCREENFITS('Property','Value',...) creates a new SCREENFITS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ScreenFits_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ScreenFits_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ScreenFits

% Last Modified by GUIDE v2.5 06-Jul-2017 09:53:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ScreenFits_OpeningFcn, ...
    'gui_OutputFcn',  @ScreenFits_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ScreenFits is made visible.
function ScreenFits_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ScreenFits (see VARARGIN)

% Choose default command line output for ScreenFits
handles.output = hObject;

% store the default cursor for interacting with the images
originalPointer = getptr(hObject);
setappdata(hObject,'originalPointer',originalPointer);% Update handles structure

% 
% % draw the big sphere in the background for rotating cells
% handles.bigAxes = axes('parent',handles.uipanel1,'units','normalize','position',[0,0,1,1]);
% [x,y,z] = sphere(15);
% % sphrScl = str2double(get(handles.edit4,'String'));
% % use a hard coded size of the sphere, not sure why it would need to be
% % something different? BPB 20151105
% sphrScl = 1e3;
% 
% % display the sphere
% surf(sphrScl*x,sphrScl*y,sphrScl*z,'Parent',handles.bigAxes,'EdgeColor','none');
% alpha(get(handles.bigAxes,'Children'),0);
% % turn it into the proper aspect ratio
% axis(handles.bigAxes,'equal');
% % set up the rotation link
% rotateLink = linkprop(handles.bigAxes,...
%     {'CameraPosition','CameraUpVector','CameraTarget'});
% setappdata(handles.figure1,'rotateLink',rotateLink);
% zoom(handles.bigAxes,4);
handles = newBigAxesWithSphere(handles);

% assign new functions for rotate3d mode
handles.rotateTool = rotate3d(handles.figure1);
% set(handles.rotateTool,'ActionPreCallback',{@hideButtonsDuringRotate,handles});
set(handles.rotateTool,'ActionPostCallback',{@postRotateFcn,handles});


% check to see if we are in HG1 or HG2
isHG1 = strfind(get(gcf,'javaFrame'),'hg.peer.HG1');
if isempty(isHG1)
    isHG1 = 0;
else
    isHG1 = 1;
end

listOfObjects = fieldnames(handles);
% for ii = 1:length(listOfObjects)
%     set(handles,listOfObjects(ii),'Units','Pixels');
% end

if isHG1
   set(handles.checkbox3,'ForeGroundColor',[0.5,0.5,0.5],'value',1); 
else
    set(handles.checkbox3,'ForeGroundColor','black','value',0);
end
guidata(hObject, handles);

% UIWAIT makes ScreenFits wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = ScreenFits_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
currentImageIdx = getappdata(handles.figure1,'imageStartIdx');
setappdata(handles.figure1,'currentImageIdx',currentImageIdx);

displayShapes(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

currentImageIdx = getappdata(handles.figure1,'imageStartIdx');
setappdata(handles.figure1,'currentImageIdx',currentImageIdx);

displayShapes(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function oldSize = numberCurrentFigs(handles)
% loop over each of elements of child handles, parse their name and make
% sure all are present
childHandles = get(handles.uipanel1,'Children');
cellHandles = cell2mat(struct2cell(handles));
storedFields = fieldnames(handles);

try
    
    for iiChild = 1:numel(childHandles);
        % strangely, the tags seem to be lost, one has to find the correct
        % name another way
        childInd = cellHandles==childHandles(iiChild);
        childName = storedFields(childInd);
        
        % hard coded split location, assuming name is axesXXX_YYY
        rowInd = str2double(childName(5:7));
        colInd = str2double(childName(9:11));
        
        % add a 1 the proper location
        axesExistArr(rowInd,colInd) = true;
        
    end
    
    if any(axesExistArr(:))
        oldSize = [0,0]; % delete and remake
    else
        oldSize = [size(axesExistArr)]; % assumes 2D array
    end
catch ME
    oldSize = [0,0]; % delete and remake
end


% --- Executes on button press in displayShapes.
function displayShapes(hObject, eventdata, handles)
% hObject    handle to displayShapes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% display the current status
set(handles.text6,'string',...
    ['File ',num2str(getappdata(handles.figure1,'currentImageIdx')),' of ',...
    num2str(getappdata(handles.figure1,'totalFrames')),'.']);

nRows = str2double(get(handles.edit1,'string'));
nCols = str2double(get(handles.edit2,'string'));

% get a list of all the children of the panel object
childHandles = get(handles.uipanel1,'Children');
boxHandles = getappdata(handles.figure1,'boxHandles');
axesHandles = getappdata(handles.figure1,'axesHandles');
axesHandles_names = getappdata(handles.figure1,'axesHandles_names');
boxHandles_names = getappdata(handles.figure1,'boxHandles_names');
rotateLink = getappdata(handles.figure1,'rotateLink');

% check to see if they are arrayed correctly (reusing axes is faster than
% recreating them

% default is to make new, then, if they happen to be correct already, we
% can ignore it

try
    for iiChild = numel(axesHandles):-1:1;
        % parse the name of this axes object
        %         % get a list of all of the fields, not by name but the actual
        %         % handle objects
        %         cellHandles = struct2cell(handles);
        %         % get a list of all of the fields, this is all the HG objects
        %         storedFields = fieldnames(handles);
        %         %
        %         childInd = cellHandles==childHandles(iiChild);
        childName = axesHandles_names{iiChild};
        % need to check if this name refers to an hg object that still
        % exists
        if ishghandle(axesHandles(iiChild));
            
            % hard coded split location to parse out row and column, assuming name is axesXXX_YYY
            rowInd = str2double(childName(5:7));
            colInd = str2double(childName(9:11));
            
            % add a 1 the proper location
            axesExistArr(rowInd,colInd) = true;
        else % no longer a valid handle
            axesHandles(iiChild) = [];
            boxHandles(iiChild) = [];
        end
    end
    
    if any(not(axesExistArr(:)))
        oldSize = [0,0]; % delete and remake
    else
        oldSize = [size(axesExistArr)]; % assumes 2D array
    end
catch ME
    oldSize = [0,0]; % delete and remake
end

if or(oldSize(1)~=nRows,oldSize(2)~=nCols)
    isNeedsNewAxes = true;
else
    isNeedsNewAxes = false;
end

if isNeedsNewAxes
    % remove any existing links to axes in the rotation
    if not(isempty(axesHandles))
        for jjAxes = 1:length(axesHandles)
            removetarget(rotateLink,axesHandles(jjAxes));
        end
    end
    % clear variables to be ready for new values
    clear('axesHandles','boxHandles','axesHandles_names','boxesHandles_names');

    % delete previous children to get ready for new ones
    for iiChild = 1:numel(childHandles)
        % remove their fields from handles
        if isfield(handles,get(childHandles(iiChild),'Tag'))
            rmfield(handles,get(childHandles(iiChild),'Tag'));
        end
        
        % delete the actual axes
        delete(childHandles(iiChild));
    end
    
    % axes height and width
    % 5% of each dimension is separation between axes
    
    axWidth = (1-0.05*(nCols+1) )/(nCols);
    axHeight = (1-0.05*(nRows+1))/(nRows);
    counter=1;
    
    for jjRow = 1:nRows
        for kkCol = 1:nCols
            axName = ['axes',sprintf('%03d',jjRow),'_',sprintf('%03d',kkCol)];
            boxName=['box',sprintf('%03d',jjRow),'_',sprintf('%03d',kkCol)];
            handles.(axName) = axes('parent',handles.uipanel1,...
                'Units','normalized',...
                'Position',...
                [axWidth*(kkCol-1)+0.05*kkCol,...
                axHeight*(jjRow-1)+0.05*jjRow,...
                axWidth,...
                axHeight]);
            
                addtarget(rotateLink,handles.(axName));
            
                % make a different box depending on HG1 or HG2
                if get(handles.checkbox3,'Value')
                 handles.(boxName)=uicontrol('parent',handles.uipanel1,...
                'style','checkbox','units','Normalized',...
                'position',[axWidth*(kkCol-1)+0.05*kkCol,...
                axHeight*(jjRow-0.9)+0.05*jjRow,...
                min(axWidth,0.1),...
                min(axHeight,0.04)],'string','flag','visible','on');
                else
            handles.(boxName)=uicontrol('parent',handles.uipanel1,...
                'style','checkbox','units','Normalized',...
                'position',[axWidth*(kkCol-1)+0.05*kkCol,...
                axHeight*(jjRow-1)+0.05*jjRow,...
                axWidth,...
                axHeight],'string','flag','visible','on');
                % find the underlying java swing object and set its opacity to
            % off
            jSwingHand = findjobj(handles.(boxName));
            setappdata(handles.(boxName),'jSwingHand',jSwingHand);
            jSwingHand.setOpaque(0);
            
                end
            
            % build the box and the axes, but don't plot anything yet
            axesHandles_names{counter}=axName;
            axesHandles(counter) = handles.(axName);
            boxHandles(counter)=handles.(boxName);
            boxHandles_names{counter} = boxName;
            counter=counter+1;
            set(handles.(axName),...
                'Tag',axName);
        end
        
    end
    % check to see if the big sphere is still there
    handles = newBigAxesWithSphere(handles);
else % reset all of the checkboxes back to false
    for jjRow = 1:nRows
        for kkCol = 1:nCols
            boxName=['box',sprintf('%03d',jjRow),'_',sprintf('%03d',kkCol)];
            set(handles.(boxName),'Value',get(handles.(boxName),'Min'));
        end
    end
    
end

for jjRow = 1:nRows
    for kkCol = 1:nCols
        axName = ['axes',sprintf('%03d',jjRow),'_',sprintf('%03d',kkCol)];
        if  getappdata(handles.figure1,'currentImageIdx')<=getappdata(handles.figure1,'totalFrames');
            plotter2(handles.(axName));
        else
            cla(handles.(axName));
            boxName=['box',sprintf('%03d',jjRow),'_',sprintf('%03d',kkCol)];
            set(handles.(boxName),'Value',get(handles.(boxName),'Min'));
        end
    end
end
colormap(handles.figure1,flipud(cbrewer('div','PiYG',256)));
set(axesHandles,'clim',[-15e-6,15e-6]);
setappdata(handles.figure1,'axesHandles',axesHandles);
setappdata(handles.figure1,'boxHandles',boxHandles);
setappdata(handles.figure1,'boxHandles_names',boxHandles_names);
setappdata(handles.figure1,'axesHandles_names',axesHandles_names);
setappdata(handles.figure1,'rotateLink',rotateLink);

% grab the children, now they should be the correct ones
childHandles = get(handles.uipanel1,'Children');
guidata(handles.figure1,handles);



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)

switch eventdata.Key
    case 'alt'
        if strcmp(get(handles.uitoggletool1,'state'),'off')
            set(handles.uitoggletool1,'state','on');
        else
            set(handles.uitoggletool1,'state','off')
        end
%         setptr(hObject,'rotate');
%         set(hObject,'WindowKeyPressFcn',[]);
%         set(handles.rotateTool,'enable','on');
    case 'leftarrow'
%         if strcmp(eventdata.Modifier,'command')
            % only use left/right when not in rotate mode
            if strcmp(get(handles.uitoggletool1,'state'),'off')
            prevImages_Callback(hObject, eventdata, handles)
            end
%         end
    case {'rightarrow'}
       % only use left/right when not in rotate mode
       if strcmp(get(handles.uitoggletool1,'state'),'off')
%         if strcmp(eventdata.Modifier,'command')
            nextImages_Callback(hObject, eventdata, handles)
        end
    otherwise
%         disp(eventdata.Character)
%         disp(eventdata.Key)
        
end



% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in selectFolder.
function selectFolder_Callback(hObject, eventdata, handles)
% hObject    handle to selectFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% hObject    handle to pbSelectFolders (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% persist startup location to base matlab
if isappdata(0,'mostRecentFolder');
    startLocation = getappdata(0,'mostRecentFolder');
else
    startLocation = [];
end
setappdata(handles.figure1,'fileEnder',[get(handles.fileEnder,'string') '.mat']);

% replace call to uipickfiles with uigetdir
folderlist = uigetdir(startLocation);
% check to make sure there was a folder selected
if not(folderlist)
    return    
end
folderlist = {folderlist};
setappdata(handles.figure1,'folder',folderlist{1});
setappdata(handles.figure1,'currentFolder',1);
set(handles.text5,'string',folderlist{1});
%set(handles.folderName,'String',folderlist{1});
% most recent folder is the parent of the last file in the list
[mostRecentFolder,~] = fileparts(folderlist{end});
% [mostRecentFolder,~] = fileparts(mostRecentFolder); % grandparent
setappdata(0,'mostRecentFolder',mostRecentFolder);

d=dir([folderlist{1} filesep '*' get(handles.fileEnder,'string') '.mat']);
if isempty(d)
    warning('ScreenFits:noFiles','No files selected')
    return
end
setappdata(handles.figure1,'imageList',d);
setappdata(handles.figure1,'totalFrames',length(d));
setappdata(handles.figure1,'currentImage',d(1).name);
setappdata(handles.figure1,'currentImageIdx',1);
setappdata(handles.figure1,'imageStartIdx',1);
set(handles.text6,'string',['File 1 of', num2str(length(d)),'.']);

displayShapes(hObject, eventdata, handles);



function fileEnder_Callback(hObject, eventdata, handles)
% hObject    handle to fileEnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileEnder as text
%        str2double(get(hObject,'String')) returns contents of fileEnder as a double


% --- Executes during object creation, after setting all properties.
function fileEnder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileEnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nextImages.
function nextImages_Callback(hObject, eventdata, handles)
% hObject    handle to nextImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check to see if this is past the end
if getappdata(handles.figure1,'currentImageIdx')==-1
    return
end

% flag the cells that have been selected
nRemoved = flagCells(handles);

% update the index for which cells to display. One can't just increment by
% a whole screen's worth because there may be files missed once some are
% removed.
nImages  = str2double(get(handles.edit1,'String'))*str2double(get(handles.edit2,'String'));
currentImageIdx = getappdata(handles.figure1,'currentImageIdx');

currentImageIdx = currentImageIdx-nRemoved;
setappdata(handles.figure1,'currentImageIdx',currentImageIdx);
setappdata(handles.figure1,'imageStartIdx',currentImageIdx);
axesHandles = getappdata(handles.figure1,'axesHandles');
boxHandles = getappdata(handles.figure1,'boxHandles');

% if we are at the end of all the files, clear all the axes so that the
% user knows she is finished
if    getappdata(handles.figure1,'currentImageIdx')>getappdata(handles.figure1,'totalFrames')
    set(handles.text6,'string',['Finished all ',num2str(getappdata(handles.figure1,'totalFrames')), ' files.']);
    for iPanel=1:length(axesHandles)
        set(boxHandles(iPanel),'Value',get(boxHandles(iPanel),'Min'));
        cla(axesHandles(iPanel));
        
    end
    setappdata(handles.figure1,'currentImageIdx',-1);
else
    % move on to the next set
    displayShapes(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function uitoggletool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentState = get(handles.uitoggletool1,'state');
switch currentState
    case 'off'
        rotate3d(handles.rotateTool,'off');
    case 'on'
        rotate3d(handles.rotateTool,'on');
end

% % % handles.rotateTool = rotate3d(handles.figure1);
% % % guidata(handles.figure1,handles);
% % % % keyboard;
% % % % look to see if there are any boxes
% % % 
% % % if any(strncmp(fieldnames(handles),'box',3))
% % %     % pull out the box handles
% % %     boxHandles_names =  getappdata(handles.figure1,'boxHandles_names');
% % %     boxHandles = getappdata(handles.figure1,'boxHandles');
% % % end
% % % 
% % % % check to see if one needs to draw the big sphere in the background
% % % drawNew = not(isfield(handles,'bigAxes'));
% % % if not(drawNew)
% % %     drawNew = not(ishghandle(handles.bigAxes));
% % % end
% % % 
% % % if drawNew
% % %     handles.bigAxes = axes('parent',handles.uipanel1,'units','normalize','position',[0,0,1,1]);
% % %     [x,y,z] = sphere(15);
% % %     sphrScl = str2double(get(handles.edit4,'String'));
% % %     %    sphrScl = 1e3;
% % %     surf(sphrScl*x,sphrScl*y,sphrScl*z,'Parent',handles.bigAxes,'EdgeColor','none');
% % %     
% % %     axis(handles.bigAxes,'equal');
% % %     rotateLink = getappdata(handles.figure1,'rotateLink');
% % %     addtarget(rotateLink,handles.bigAxes);
% % %     
% % %     zoom(handles.bigAxes,4);
% % %     guidata(handles.figure1,handles);
% % %     
% % % end
% % % 
% % % if strcmp(get(handles.rotateTool,'enable'),'on');
% % % %     setappdata(handles.figure1,'isRotateMode',false)
% % %     set(boxHandles,'visible','on');
% % %     set(boxHandles,'enable','on');
% % %     set(handles.bigAxes,'visible','off');
% % %     %     set(get(handles.bigAxes,'Children'),'visible','off');
% % %     set(handles.rotateTool,'enable','off');
% % %     % prevent accidental button pushing
% % %     set(handles.selectFolder,'enable','on','visible','on');
% % %     set(handles.prevImages,'enable','on','visible','on');
% % %     set(handles.nextImages,'enable','on','visible','on');
% % %     set(handles.popupmenu1,'enable','on','visible','on');
% % %     
% % % elseif strcmp(get(handles.rotateTool,'enable'),'off');
% % %     set(boxHandles,'enable','off');
% % %     set(boxHandles,'visible','off');
% % %     %     set(get(handles.bigAxes,'Children'),'visible','on');
% % %     axis(handles.bigAxes,'off');
% % %     %     set(handles.bigAxes,'grid','off');
% % %     alpha(get(handles.bigAxes,'Children'),0);
% % %     set(handles.rotateTool,'enable','on');
% % %     % prevent accidental button pushing
% % %     set(handles.selectFolder,'enable','off','visible','off');
% % %     set(handles.nextImages,'enable','off','visible','off');
% % %     set(handles.prevImages,'enable','off','visible','off');
% % %     set(handles.popupmenu1,'enable','off','visible','off');
% % %     setappdata(handles.figure1,'isRotateMode',false);
% % % 
% % % else
% % %     %     disp('something else');
% % % end

function handles = newBigAxesWithSphere(handles)
% check to see if one needs to draw the big sphere in the background
drawNew = not(isfield(handles,'bigAxes'));
if not(drawNew)
    drawNew = not(ishghandle(handles.bigAxes));
%     if isappdata(handles.figure1,'rotateLink')
%         rotateLink = getappdata(handles.figure1,'rotateLink');
%         rmtarget(rotateLink, handles.bigAxes);
%     end
end

if drawNew
    handles.bigAxes = axes('parent',handles.uipanel1,'units','normalize','position',[0,0,1,1]);
    [x,y,z] = sphere(15);
    sphrScl = 1e3;
    surf(sphrScl*x,sphrScl*y,sphrScl*z,'Parent',handles.bigAxes,'EdgeColor','none');
    alpha(get(handles.bigAxes,'Children'),0);
    set(handles.bigAxes,'Visible','off');
    axis(handles.bigAxes,'equal');
    if isappdata(handles.figure1,'rotateLink');
     rotateLink = getappdata(handles.figure1,'rotateLink');
         addtarget(rotateLink,handles.bigAxes);
    else
        rotateLink = linkprop(handles.bigAxes,...
            {'CameraPosition','CameraUpVector','CameraTarget'});
        
    end
    try
    zoom(handles.bigAxes,4);
    catch ME % sometimes (specifically R2016a) this throws a reshape error
        zoom(handles.bigAxes,4);
    end
    guidata(handles.figure1,handles);
end

setappdata(handles.figure1,'rotateLink',rotateLink);

function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in deleteImages.
function deleteImages_Callback(hObject, eventdata, handles)
% hObject    handle to deleteImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of deleteImages


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% currentImageIdx = getappdata(handles.figure1,'currentImageIdx');
% nImages  = str2double(get(handles.edit1,'String'))*str2double(get(handles.edit2,'String'));
% currentImageIdx = max([1,currentImageIdx-1*nImages]);
currentImageIdx = getappdata(handles.figure1,'imageStartIdx');
setappdata(handles.figure1,'currentImageIdx',currentImageIdx);

displayShapes(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(hObject,'Value') == get(hObject,'Max');
    h = msgbox('Enabling ''drawnow'' will update the display at the end of each file but may decrease the overall performance.');
end


% --- Executes on button press in prevImages.
function prevImages_Callback(hObject, eventdata, handles)
% hObject    handle to prevImages (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axesHandles=getappdata(handles.figure1,'axesHandles');
boxHandles=getappdata(handles.figure1,'boxHandles');
axesHandles_names=getappdata(handles.figure1,'axesHandles_names');
boxHandles_names=getappdata(handles.figure1,'boxHandles_names');

% check to see if this is past the end
if getappdata(handles.figure1,'currentImageIdx')==-1
    
    nImages  = str2double(get(handles.edit1,'String'))*str2double(get(handles.edit2,'String'));
    totalFrames = getappdata(handles.figure1,'totalFrames');
    setappdata(handles.figure1,'currentImageIdx',max(1,totalFrames-nImages+1))
    
    displayShapes(hObject, eventdata, handles);
    return
end

% check to see if there were any cells selected before going back a page
isABoxChecked = false;
iPanel = 1;
while not(isABoxChecked)
    isABoxChecked = get(boxHandles(iPanel),'Value');
    if iPanel == length(boxHandles)
        break
    end
    iPanel = iPanel+1;
end

% if there was a cell selected, prompt the user for what they would like to
% do
if isABoxChecked
    answer1 = questdlg('Do you want to flag selected cells?','Save flags','Yes and go back a page','No, but still go back a page','Cancel','Yes and go back a page');
else
    answer1 = 'No, but still go back a page';
end

% process their decision to save and go back, go back or cancel navigation
switch answer1
    case 'Yes and go back a page'
        % apply the flag to selected images
        nRemoved = flagCells(handles);
        
        % move back to the previous page
        currentImageIdx = getappdata(handles.figure1,'currentImageIdx');
        % calculate the number of images displayed on one "page"
        nImages  = str2double(get(handles.edit1,'String'))*str2double(get(handles.edit2,'String'));
        % move the image index back one page
        currentImageIdx = max([1,currentImageIdx-2*nImages]);
        
        setappdata(handles.figure1,'currentImageIdx',currentImageIdx);
        setappdata(handles.figure1,'imageStartIdx',currentImageIdx);

    case 'No, but still go back a page'
        % move back to the previous page
        currentImageIdx = getappdata(handles.figure1,'currentImageIdx');
        % calculate the number of images displayed on one "page"
        nImages  = str2double(get(handles.edit1,'String'))*str2double(get(handles.edit2,'String'));
        % move the image index back one page
        currentImageIdx = max([1,currentImageIdx-2*nImages]);
        setappdata(handles.figure1,'currentImageIdx',currentImageIdx);
        setappdata(handles.figure1,'imageStartIdx',currentImageIdx);

    otherwise
        
        return
end

displayShapes(hObject, eventdata, handles);

function nRemoved = flagCells(handles)
% this function applies the the "FLAG" appendix to all the files that are
% flagged in the current set of displayed images. If the user desires, it
% also deletes the tifs associated with these files

% this function should be called when moving forward or backward to a new
% set of files

% obtain stored data
axesHandles = getappdata(handles.figure1,'axesHandles');
boxHandles = getappdata(handles.figure1,'boxHandles');
imageList = getappdata(handles.figure1,'imageList');
currentImageIdx = getappdata(handles.figure1,'currentImageIdx');
totalFrames = getappdata(handles.figure1,'totalFrames');
nImagesPerField = length(axesHandles);

% counter for the number of files that have been removed
nRemoved = 0;
% loop over all the files that need to be moved
for i=length(axesHandles):-1:1
    if    get(boxHandles(i),'Value')
        % increment counter
        nRemoved = nRemoved+1;
        
        % retrieved the name of the current file
        fname=getappdata(axesHandles(i),'currentImage');
        
        % replace its appendix with a flag
        fname2=strrep(fname,'.mat','FLAG.mat');
        
        % move file
        movefile(fname,fname2)
        
        % remove it from the list of files
        imageList(currentImageIdx+i-nImagesPerField-1) = [];
        
        % decrement the count of total frames
        totalFrames = totalFrames - 1;
        
        % if the user wants the tifs deleted, do so
        if get(handles.deleteImages,'Value')
            tiffname=strrep(fname,[get(handles.fileEnder,'string') '.mat'],'.tif');
            delete(tiffname);
        end
    end
end

% store the modified list of files and counter for image index
setappdata(handles.figure1,'imageList',imageList)
setappdata(handles.figure1,'totalFrames',totalFrames);
set(handles.text6,'string',...
    ['File ',num2str(getappdata(handles.figure1,'currentImageIdx')),' of ', ...
    num2str(getappdata(handles.figure1,'totalFrames')),'.']);


% --- Executes on key release with focus on figure1 or any of its controls.
function figure1_WindowKeyReleaseFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was released, in lower case
%	Character: character interpretation of the key(s) that was released
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) released
% handles    structure with handles and user data (see GUIDATA)

function figure1_KeyReleaseFcn(hObject, eventdata, handles)
% this function is apparently expected, but it shouldn't do anything


% --------------------------------------------------------------------
function uitoggletool1_OffCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
boxHandles = getappdata(handles.figure1,'boxHandles');
originalPointer = getappdata(handles.figure1,'originalPointer');
set(handles.figure1,originalPointer{:});

% display the boxes for check marks and enable them
set(boxHandles,'visible','on');
set(boxHandles,'enable','on');
set(handles.bigAxes,'visible','off');
   
% renable button pushing
set(handles.selectFolder,'enable','on','visible','on');
set(handles.prevImages,'enable','on','visible','on');
set(handles.nextImages,'enable','on','visible','on');
set(handles.popupmenu1,'enable','on','visible','on');

rotate3d(handles.figure1,'off');
% hManager = uigetmodemanager(handles.figure1);
% try
%     set(hManager.WindowListenerHandles, 'Enable', 'on');  % HG1
% catch
%     [hManager.WindowListenerHandles.Enabled] = deal(true);  % HG2
% end



% --------------------------------------------------------------------
function uitoggletool1_OnCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

boxHandles = getappdata(handles.figure1,'boxHandles');
% disable view of the checkboxes and display invisbly the rotate axis
set(boxHandles,'enable','off');
set(boxHandles,'visible','off');
axis(handles.bigAxes,'off');
alpha(get(handles.bigAxes,'Children'),0);

% prevent accidental button pushing
set(handles.selectFolder,'enable','off','visible','off');
set(handles.nextImages,'enable','off','visible','off');
set(handles.prevImages,'enable','off','visible','off');
set(handles.popupmenu1,'enable','off','visible','off');
set(handles.rotateTool,'enable','on');

function postRotateFcn(hObj,eventdata,handles)

hManager = uigetmodemanager(handles.figure1);
try
    set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
catch
    [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
end

set(handles.figure1, 'WindowKeyPressFcn', {@figure1_WindowKeyPressFcn,handles});


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4



axesHandles = getappdata(handles.figure1,'axesHandles');
boxHandles = getappdata(handles.figure1,'boxHandles');

if get(hObject,'Value')
   axis(axesHandles,'vis3d'); 
else
    axis(axesHandles,'equal');
end
