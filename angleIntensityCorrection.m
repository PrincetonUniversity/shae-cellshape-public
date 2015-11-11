function [thetaCorrection,thetaAngleOut]=angleIntensityCorrection

% this is from WT mreB cells
%    thetaI= [0.2908    0.2916    0.2966    0.3026    0.3113    0.3221    0.3349 ...
%        0.3457    0.3525    0.3628    0.3767 0.3803    0.3769    0.3718    0.3687...
%        0.3553    0.3432    0.3307    0.3210    0.3099    0.3000    0.2965 ...
%        0.2908    0.2892    0.2892    0.2908    0.2965    0.3000    0.3099 ...
%        0.3210    0.3307    0.3432    0.3553 0.3687    0.3718    0.3769  ...
%        0.3803    0.3767    0.3628    0.3525    0.3457    0.3349    0.3221  ...
%        0.3113 0.3026    0.2966    0.2916    0.2908];
%    
%    thetaAngleOut=linspace(-pi/2,3*pi/2,length(thetaI));
%    
   
% this is from 1ug/mL mreB cells
 thetaI=   [0.2897    0.2941    0.2858    0.3044    0.3117    0.3149    0.3371    0.3502    0.3472    0.3592    0.3800 ...
0.3871    0.3687    0.3678    0.3717    0.3575    0.3241    0.3350    0.3200    0.2972    0.3036    0.2914 ...
0.2916    0.2899];

thetaI=[fliplr(thetaI),thetaI,fliplr(thetaI)];
thetaI=smooth(thetaI,5);
thetaI=thetaI(length(thetaI)/3:end);
   thetaAngleOut=linspace(-pi/2,pi/2,length(thetaI));


   
   
thetaI=thetaI/mean(thetaI);
thetaCorrection=1./thetaI';