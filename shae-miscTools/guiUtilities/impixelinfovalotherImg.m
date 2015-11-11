
% % BPB 2012-10-13 calls impixelinfoval to reference a different image

function huic = impixelinfovalotherImg(parent,himage,otherImg)

% use the builtin impixelinfo to find the position of the cursor but
% display the information from a different image

huic = impixelinfoval(parent,himage);
iptaddcallback(parent,'WindowButtonMotionFcn',{@updateTextString,huic,otherImg});
% javaHandle = findjobj(huic);
% set(javaHandle, 'PropertyChangeCallback',{@updateTextString,huic,otherImg});


function updateTextString(hobj,evnt,huic,otherImg)
% findhandles
handles = guidata(get(huic,'parent'));

originalString = get(huic,'String');

newString = originalString;
% this may have issues if the default text in impixelvalinfo chages
if strncmp(originalString,'(X, Y)',6)
   newString(9:end) = [];
   newString = [newString,'[NaN]'];
else
   % find first comma  separating the x and y coordinates
   indComma = strfind(originalString,',');
   colInd = str2double(originalString(2:indComma(1)-1));
   
   % find first ], the end of the y coordinate
   indParen = strfind(originalString,')');
   rowInd = str2double(originalString(indComma(1)+2:indParen(1)-1));
  
    % pull out the value of the otherImage at this point
    imgMod2 = imagemodel(otherImg);
%     pixelValString = getPixelInfoString(imgMod2,rowInd,colInd);
    pixelVal = getPixelValue(imgMod2,rowInd,colInd);
    newString(indParen+1:end) = [];
    newString = [newString,' [',num2str(pixelVal),']'];
    
end

set(huic,'String',newString);

