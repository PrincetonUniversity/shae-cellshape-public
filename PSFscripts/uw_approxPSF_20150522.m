% membrane width in images appears to be about 4 pixels fwhm
psfWidthXY = 1.7; 
% membrane width in images appears to be about 2.5 pixels fwhm (125 nm
% steps)
psfWidthZ = 1;

psfXY = exp(-([-20:20]).^2./(2*psfWidthXY.^2));
psfZ =  exp(-([-40:40]).^2./(2*psfWidthZ.^2));

[psfX,psfY,psfZ] = ndgrid(psfXY,psfXY,psfZ);

psfApprox = psfX.*psfY.*psfZ;
psfApprox = normalizeRange(psfApprox);
 setappdata(0,'averageSurfPSF',psfApprox);
