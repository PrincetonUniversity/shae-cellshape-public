function [image_para,flag,cline_para] = paramsText(fullParamsFile);

[cline_para,image_para,flag] = parseParams(fullParamsFile);

% pull out variables that are used outside of structures
% ----
pixel_interp = flag.pixel_interp;
intensity_correction = flag.intensity_correction;
window = flag.window;
end_angles = flag.end_angles;
pixel_interp = flag.pixel_interp;
Fflag2 = flag.Fflag2;

% ----

%% manipulate primary parameters into secondary parameters
try
% image parameters
image_para.ZF_factor=...
    image_para.Z_scale*image_para.Fstack_seperation_nm/image_para.nm_per_pixel;
image_para.nm_per_zstack=...
image_para.Z_scale*image_para.stack_seperation_nm;
image_para.ZX_factor=image_para.nm_per_zstack/image_para.nm_per_pixel;

% cline paramaters
if rem(cline_para.pts,4)~=0;
    % warn about improper size
    warning('parseParams:cline_para_ptsMod4',...
        'Changing cline_para.pts to be divisble by 4');
    cline_para.pts=cline_para.pts+4-rem(cline_para.pts,4);
end

cline_para.radius=cline_para.radius*pixel_interp; %initial radius
cline_para.center=[pixel_interp*window+2,pixel_interp*window+2];

%% % ----------- begin default params
if ~isfield(cline_para,'alpha') && isfield(cline_para,'stiff')
    cline_para.alpha =1*cline_para.stiff; 
    cline_para.beta = cline_para.pts*cline_para.stiff;
end

if ~isfield(cline_para,'alpha3d') && isfield(cline_para,'stiff3d')
    cline_para.alpha3d =1*cline_para.stiff3d; 
    cline_para.beta3d = cline_para.pts*cline_para.stiff3d; 
end

if ~isfield(cline_para,'zalpha3d')
    cline_para.zalpha3d=cline_para.alpha3d;
    cline_para.zbeta3d=cline_para.beta3d;
end    

if ~isfield(cline_para,'gamma')
cline_para.gamma=10/pixel_interp; %step size
end

% settings for interpolating fluorescence after the fact
if ~isfield(flag,'imageOutlineStack')
    flag.imageOutlineStack = 0;
end

if ~isfield(flag,'shapeName')
    flag.shapeName = 'c1';
end

if ~isfield(flag,'proteinName')
    flag.proteinName = 'c2';
end

if ~isfield(flag,'profile')
    flag.profile = '0';
end

% % % ----------- end default params
%%
% check ordering of stack
if flag.reverse==1 && flag.F1==0
    warning('parseParams:reverseNoF1',...
        ['No protein channel selected but shape channel is set as last channel. ',...
        'Assuming one wasted channel and then a shape channel.']);
flag.F1 = 1;
end



eval(flag.psfScript);

catch ME % go ahead and store the workspace, make sure to clean this file up later
   save('errorSpace.mat');
   rethrow(ME);
end



%% these appear to be unused
% t_start=1;          %start and end time
% t_end=1;
% if t_end>image_para.stack_t_size
%     t_end=image_para.stack_t_size;
% end

function  [cline_para,image_para,flag] = parseParams(textFile);
%%
fidA = fopen(textFile);

while not(feof(fidA))
    nextLine = fgetl(fidA);
    % is the line blank? skip it
    if isempty(nextLine)
        continue
    end
    
    % is the first character a '%', skip it
    if strcmp(textscan(nextLine,'%c',1),'%')
        continue
    end
    
    % is the first word the beginning of a section?
    [firstWord,strpos] = textscan(nextLine,'%s',1);
    firstWord = firstWord{1}{1};
    nextLine = nextLine(strpos:end);
    if strcmp(firstWord,'beginSection')
        % pull out the name of the struct
        structName = textscan(nextLine,'%s',1);
        structName = structName{1}{1};
        tempStruct = struct();
        continue
    end
    
    if strcmp(firstWord,'endSection')
        % store the structure
        eval([structName,' = tempStruct;']);
        continue
    end
    variableName = firstWord;
    [value,strpos] = textscan(nextLine,'%f',1);
    value = value{1};
    % flag for string values;
    if value==-999
        nextLine = nextLine(strpos:end);
        [value,strpos] = textscan(nextLine,'%s',1);
        value = value{1}{1};
    else
    end
    
    % store the value
    tempStruct.(variableName) = value;
end

% warn if variables don't exist
if not(exist('cline_para','var'));
    warning('parseParams:cline_para','Variable ''cline_para'' does not exist');
    cline_para = struct([]);
end
if not(exist('image_para','var'));
    warning('parseParams:image_para','Variable ''image_para'' does not exist');
    image_para = struct([]);
end
if not(exist('flag','var'));
    warning('parseParams:flag','Variable ''flag_para'' does not exist');
    image_para = struct([]);
end

fclose(fidA);


