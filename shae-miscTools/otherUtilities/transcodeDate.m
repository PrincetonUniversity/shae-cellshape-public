function a = transcodeDate(b)
% This function has two types, decode and encode

%% inputs 
% b is either the string or serial number to be coded or decoded

%% outputs
% a is of the other format

%% Code
lut = ['0','1','2','3','4','5','6','7','8','9',...
        'a','b','c','d','e','f','g','h','i','j',...
        'k','l','m','n','o','p','q','r','s','t',...
        'u','v','w','x','y','z','A','B','C','D',...
        'E','F','G','H','I','J','K','L','M','N',...
        'O','P','Q','R','S','T','U','V','W','X','Y','Z'];

if isnumeric(b) % encode serial number to string
    timeStr = '00000000';
    timeSN = b*10^10;
for iChar = 1:9;
    timeStr(iChar) = lut(1+floor(timeSN/62^(9-iChar)));
    timeSN = rem(timeSN,62^(9-iChar));
end
    a = timeStr;
    return
else
    if ischar(b) % characters to decode
        if and(eq(size(b,2),20),strfind(b,':')) % decoding andor timestamp
            a = datenum(b,'yyyy:mm:dd HH:MM:SS');
          return
        else % decode encoded string to serial number
            
%             datestr
            timeSN = 0;
            timeStr = b;
            
            for iChar = 9:-1:1
                timeSN = timeSN + (-1+findstr(timeStr(iChar),lut))*(62^(9-iChar));
            end
            timeSN = timeSN/(10^10);
            a = timeSN;
            return
        end
        
        else
            a = false;
            return
        end
    end
end