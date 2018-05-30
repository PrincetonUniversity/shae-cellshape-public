function varargout = FitViewer(varargin)
% FITVIEWER MATLAB code for FitViewer.fig
%      FITVIEWER, by itself, creates a new FITVIEWER or raises the existing
%      singleton*.
%
%      H = FITVIEWER returns the handle to a new FITVIEWER or the handle to
%      the existing singleton*.
%
%      FITVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FITVIEWER.M with the given input arguments.
%
%      FITVIEWER('Property','Value',...) creates a new FITVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FitViewer_OpeningFcn gets xrcalled.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FitViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FitViewer

% Last Modified by GUIDE v2.5 26-Feb-2014 14:52:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FitViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @FitViewer_OutputFcn, ...
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


% --- Executes just before FitViewer is made visible.
function FitViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FitViewer (see VARARGIN)
playt.TimerFcn = {@TmrFcn,handles};
playt.BusyMode = 'Queue';
playt.ExecutionMode = 'FixedRate';
playt.Period = 1/str2double(get(handles.tFPS,'String'));

setappdata(handles.figure1,'playTimer',timer(playt));

% Choose default command line output for FitViewer
handles.output = hObject;
hlistener=addlistener(handles.slider1,'ContinuousValueChange',...
    @plotter);


setappdata(handles.figure1,'holdaxes',false);
setappdata(handles.slider1,'hlistener',hlistener);
setappdata(handles.figure1,'CurrentFrame',1)
setappdata(handles.figure1,'totalFrames',14);
setappdata(handles.figure1,'m_frames',1);
setappdata(handles.figure1,'cmap','Z')

set(handles.slider1,'SliderStep',[1,1]);
set(handles.slider1,'Max',getappdata(handles.figure1,'totalFrames'));
set(handles.slider1,'Value',getappdata(handles.figure1,'m_frames'));
set(handles.slider1,'Min',getappdata(handles.figure1,'m_frames'));

figTemp = figure('visible','off');
colorMapList = {'Jet','HSV','Hot','Cool','Spring','Summer','Autumn','Winter',...
    'Bone','Copper','Pink','Parula','Black2Green'};
    
% list a set of colormaps
for ii = length(colorMapList):-1:1
    colorMapTest = colorMapList{ii};
    try
        colormap(figTemp,colorMapTest);
    catch ME
        colorMapList(ii) = [];
    end
    
end
close(figTemp)
set(handles.Cstyle,'string',colorMapList);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FitViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FitViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
%set(handles.figure1,'CurrentFrame',get(hObject,'Value'))
get(hObject,'Value');
%plotter(hObject,eventdata)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% Timer Function
function TmrFcn(src,event,handles)
% pull appdata from the handles structure
 
CurrentFrame = getappdata(handles.figure1,'CurrentFrame');
set(handles.slider1,'Value',CurrentFrame)
%CurrentFrame=1;
totalFrames = getappdata(handles.figure1,'totalFrames');
loop = false;
 
% at some point, include the ability to loop
% loop = get(handles.loop,'Value');
% if the current frame is less than the total, increment frame by one
if CurrentFrame < totalFrames
    setappdata(handles.figure1,'CurrentFrame', CurrentFrame + 1);
    set(handles.slider1,'Value',CurrentFrame+1)
    % otherwise, if looping, reset to 1
    % elseif loop == get(handles.loop,'Max')
elseif loop
    setappdata(handles.figure1,'CurrentFrame', 1);
    % otherwise, stop playback
else
    set(handles.togglePlay,'Value',get(handles.togglePlay,'Min'));
    togglePlay_Callback(handles.togglePlay, event, handles);
end
plotter(handles.slider1,event)
%guidata(handles.figure1,handles);
%updateDisplay(handles);


% --- Executes on button press in togglePlay.
function togglePlay_Callback(hObject, eventdata, handles)
% hObject    handle to togglePlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hint: get(hObject,'Value') returns toggle state of playToggle
button_state = get(hObject,'Value');
% disp(button_state);
% disp(get(hObject,'Max'));
% disp(get(hObject,'Min'));
if button_state == get(hObject,'Max')
    % Toggle button is pressed, take appropriate action
    set(hObject,'String','Stop');
    set(hObject,'ForegroundColor',[1 0 0]);
    
    %     set(handles.cursortoggle,'State','off'); % having cursor on creates errs
    start(getappdata(handles.figure1,'playTimer'))
elseif button_state == get(hObject,'Min')
    % Toggle button is not pressed, take appropriate action
    set(hObject,'String','Play');
    set(hObject,'ForegroundColor',[0 1 0]);
   
    stop(getappdata(handles.figure1,'playTimer'))
    %     set(handles.cursortoggle,'State','on'); % cursor on is default!
end


% --- Executes on button press in togglePlay.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglePlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglePlay



function tFPS_Callback(hObject, eventdata, handles)
% hObject    handle to tFPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tFPS as text
%        str2double(get(hObject,'String')) returns contents of tFPS as a double
playt = getappdata(handles.figure1,'playTimer');
set(playt,'Period',1/str2double(get(handles.tFPS,'String')));



% --- Executes during object creation, after setting all properties.
function tFPS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tFPS (see GCBO)
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
p=get(handles.edit2,'String');


cd(p)
d=dir;
l=[];
for i=1:size(d,1)
if isempty(strfind(d(i,1).name,getappdata(handles.figure1, 'fileEnder')))
l=[l,i];
end
end
d(l)=[];
setappdata(handles.figure1,'totalFrames',length(d));
set(handles.edit3,'String',num2str(length(d)));
edit3_Callback(handles.edit3,eventdata,handles);
plotter(handles.slider1,eventdata)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
set(handles.slider1,'Max',str2double(get(hObject,'String')));
set(handles.slider1,'Min',1);
set(handles.slider1,'Value',str2double(get(hObject,'String')));
setappdata(handles.figure1,'totalFrames',str2double(get(hObject,'String')));



% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
setappdata(handles.figure1,'cmap',contents{get(hObject,'Value')});
cla(handles.axes1);
plotter(handles.slider1,eventdata)
% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
setappdata(handles.figure1,'holdaxes',get(hObject,'Value'));
plotter(handles.slider1,eventdata)



function fileEnder_Callback(hObject, eventdata, handles)
% hObject    handle to fileEnder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fileEnder as text
%        str2double(get(hObject,'String')) returns contents of fileEnder as a double
setappdata(handles.figure1,'fileEnder',[get(handles.fileEnder,'string') '.mat']);
%plotter(handles.slider1,eventdata);


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

% clear the axes
cla(handles.axes1);

% replace call to uipickfiles with uigetdir
folderlist = uigetdir(startLocation);
% make sure someone actually picked a folder
if not(folderlist)
    return
end
folderlist = {folderlist};

setappdata(handles.figure1,'folderlist',folderlist{1});
setappdata(handles.figure1,'currentFolder',1);
set(handles.folderName,'String',folderlist{1});
% most recent folder is the parent of the last file in the list
[mostRecentFolder,~] = fileparts(folderlist{end});
% [mostRecentFolder,~] = fileparts(mostRecentFolder); % grandparent
setappdata(0,'mostRecentFolder',mostRecentFolder);


d=dir(folderlist{1});
l=[];
for i=1:size(d,1)
    if isempty(strfind(d(i,1).name,getappdata(handles.figure1,'fileEnder')))
        l=[l,i];
    end
end
d(l)=[];
setappdata(handles.figure1,'totalFrames',length(d));
set(handles.slider1,'Max',length(d));
set(handles.edit3,'String',num2str(length(d)));
edit3_Callback(handles.edit3,eventdata,handles);
plotter(handles.slider1,eventdata)


% --- Executes on button press in cellFlag.
function cellFlag_Callback(hObject, eventdata, handles)
% hObject    handle to cellFlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Are you sure you want to flag this Cell '  , ...
 'No', ...
 'Yes');
switch choice
    case 'Yes';
        fname=strrep(get(handles.imageFile,'string'),'.mat', '');
    movefile([getappdata(handles.figure1,'folderlist') filesep ...
        get(handles.imageFile,'string')],...
        [getappdata(handles.figure1,'folderlist') filesep ...
        fname, 'FLAG.mat'])
    otherwise
end



% --- Executes during object creation, after setting all properties.
function imageFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function trim_Callback(hObject, eventdata, handles)
% hObject    handle to trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trim as text
%        str2double(get(hObject,'String')) returns contents of trim as a double
plotter(handles.slider1,eventdata);



% --- Executes during object creation, after setting all properties.
function trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viewStack.
function viewStack_Callback(hObject, eventdata, handles)
% hObject    handle to viewStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name=[get(handles.folderName,'String') filesep get(handles.imageFile,'String')];
k=strfind(name,getappdata(handles.figure1,'fileEnder'));

frames=str2num(get(handles.FramesPerStack,'String'));
stackLoad([name(1:k-1) '.tif'],frames,.6);


% --- Executes on button press in viewFit.
function viewFit_Callback(hObject, eventdata, handles)
% hObject    handle to viewFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load([get(handles.folderName,'String') filesep get(handles.imageFile,'String')],'eg')
SliceBrowser(eg);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pbRevStart.
function pbRevStart_Callback(hObject, eventdata, handles)
% hObject    handle to pbRevStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'value',1);
plotter(handles.slider1,eventdata);

% --- Executes on button press in pbRewind.
function pbRewind_Callback(hObject, eventdata, handles)
% hObject    handle to pbRewind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'value',max(1,get(handles.slider1,'value')-5));
plotter(handles.slider1,eventdata);

% --- Executes on button press in pbBack.
function pbBack_Callback(hObject, eventdata, handles)
% hObject    handle to pbBack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider1,'value',max(1,get(handles.slider1,'value')-1));
plotter(handles.slider1,eventdata);


% --- Executes on button press in pbFwd.
function pbFwd_Callback(hObject, eventdata, handles)
% hObject    handle to pbFwd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'value',min(str2double(get(handles.edit3,'string')),...
    get(handles.slider1,'value')+1));
plotter(handles.slider1,eventdata);


% --- Executes on button press in pbFF.
function pbFF_Callback(hObject, eventdata, handles)
% hObject    handle to pbFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.slider1,'value',min(str2double(get(handles.edit3,'string')),...
    get(handles.slider1,'value')+5));
plotter(handles.slider1,eventdata);

% --- Executes on button press in pbFwdEnd.
function pbFwdEnd_Callback(hObject, eventdata, handles)
% hObject    handle to pbFwdEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.slider1,'value',str2double(get(handles.edit3,'string')));
plotter(handles.slider1,eventdata);

function jumpCurrentFrame_Callback(hObject, eventdata, handles)
% hObject    handle to jumpCurrentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of jumpCurrentFrame as text
%        str2double(get(hObject,'String')) returns contents of jumpCurrentFrame as a double

set(handles.slider1,'value',str2double(get(hObject,'string')));
plotter(handles.slider1,eventdata)



% --- Executes during object creation, after setting all properties.
function jumpCurrentFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jumpCurrentFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ImageDelete.
function ImageDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ImageDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
choice = questdlg('Are you sure you want to delete this Image? '  , ...
 'No', ...
 'Yes');
switch choice
    case 'Yes';
name=[get(handles.folderName,'String') filesep get(handles.imageFile,'String')];
k=strfind(name,'conv');
ImageName=[name(1:k-1) '.tif'];
        delete(ImageName)
        
        
        
    otherwise
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Colourbar.
function Colourbar_Callback(hObject, eventdata, handles)
% hObject    handle to Colourbar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    colorbar('peer',handles.axes1);

% Hint: get(hObject,'Value') returns toggle state of Colourbar



% --- Executes on selection change in Cstyle.
function Cstyle_Callback(hObject, eventdata, handles)
% hObject    handle to Cstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject,'String'));
Cstyle=contents{get(hObject,'Value')};
setappdata(handles.figure1,'Cstyle',Cstyle);
colormap(handles.axes1,Cstyle);
% Hints: contents = cellstr(get(hObject,'String')) returns Cstyle contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Cstyle


% --- Executes during object creation, after setting all properties.
function Cstyle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Cstyle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function FramesPerStack_Callback(hObject, eventdata, handles)
% hObject    handle to FramesPerStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FramesPerStack as text
%        str2double(get(hObject,'String')) returns contents of FramesPerStack as a double


% --- Executes during object creation, after setting all properties.
function FramesPerStack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FramesPerStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PolymerPreview.
function PolymerPreview_Callback(hObject, eventdata, handles)
% hObject    handle to PolymerPreview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
name=[get(handles.folderName,'String') filesep get(handles.imageFile,'String')];
surfCellPolymer(name);


% --- Executes on button press in rowNormalize.
function rowNormalize_Callback(hObject, eventdata, handles)
% hObject    handle to rowNormalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plotter(handles.slider1,eventdata)

% Hint: get(hObject,'Value') returns toggle state of rowNormalize
