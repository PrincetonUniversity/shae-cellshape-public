function modeName = enrichmentModeString(varargin)
% this function generates unique file names depending on the types of
% curvature used and the various other flags (undersample, midflag,
% poleflag)
% modeNameShort = enrichmnentModeString(options) returns the simplified
% string of which optiones were used
% modeNameLong = enrichmnentModeString(modeNameShort) returns the human
% readable version of the short string
% version 1.0
% version 1.1 BPB 20161101 Added background subtraction options

if nargin>1
    error('enrichmentModeString:wrongNumberInputs','Incorrect number of input parameters');
end

if isstruct(varargin{1})
    options = varargin{1};
    % short string
    if options.undersample
        modeName = '1';
    else
        modeName = '0';
    end
    
    
    
    if options.midflag
        modeName = [modeName,'1'];
    else
        modeName = [modeName,'0'];
    end
    
    if options.poleflag == 1
        modeName = [modeName,'1'];
    elseif options.poleflag == -1
        modeName = [modeName,'2'];
    else
        modeName = [modeName,'0'];
    end
    
    switch options.geometricPropertyName
        case 'gc'
            modeName = [modeName,'1'];
        case 'gm'
            modeName = [modeName,'2'];
        case 'CLT'
            modeName = [modeName,'3'];
        case 'CLK'
            modeName = [modeName,'4'];
        otherwise
            modeName = [modeName,'z'];
    end
    
    switch options.fluorescentMeasurementName
        case 'Ff'
            modeName = [modeName,'1'];
        case 'Ffr'
            modeName = [modeName,'2'];
        case 'amountNotConcentration';
            modeName = [modeName,'3'];
        otherwise
            modeName = [modeName,'z'];
    end
    
    % added in version 1.1
    switch options.backgroundSubtract
        case -999
            modeName = [modeName,'0'];      
        case -998
            modeName = [modeName,'1'];          
        case -995
            modeName = [modeName,'2'];
        case -990
            modeName = [modeName,'3'];    
        otherwise
            modeName = [modeName,'4'];         
    end
elseif ischar(varargin{1})
    % process the short form of the string into a long description
   
    modeName = varargin{1};
    switch modeName(1)
        case '1'
            disp('Undersample data');
        case '0'
            disp('Do not undersample data');
    end
    
    switch modeName(2)
        case '1'
            disp('Use only the middle ''Z'' portion of coordinates');
        case '0'
            disp('Use all ''Z'' coordinates');
    end
    
    switch modeName(3)
        case '1'
            disp('Include the poles in the data');
        case '0'
            disp('Remove poles after normalization');
        case '2'
            disp('Remove poles before normalization');
    end
    
    switch modeName(4)
        case '1'
            disp('Gaussian curvature');
        case '2'
            disp('Mean curvature');
        case '3'
            disp('Centerline torsion');
        case '4'
            disp('Centerline curvature');
        case 'z'
            disp('Unknown geometric property');
        
    end
    
    
    switch modeName(5)
        
        case '1'
            disp('Ff fluorescence at surface');
        case '2'
            disp('Ffr fluorescence max radial projection');
        case '3'
            disp('Amount of fluorescence, not concentration of fluorescence');
        case 'z'
            disp('Unknown fluorescent measurment');
    end
    
    if numel(modeName)<=5
        disp('Enrichment string version <v1.1')
    end
else
    error('enrichmentModeString:wrongInputType','Incorrect input type')
end

