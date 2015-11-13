function matrixOut = smooth2b(matrixIn,Nr,Nc)
% Smooths 2D array data.  Ignores NaN's.
%
%function matrixOut = smooth2a(matrixIn,Nr,Nc)
% 
% This function smooths the data in matrixIn using a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1).  Basically, you end up replacing
% element "i" by the mean of the rectange centered on "i".  Any NaN
% elements are ignored in the averaging.  If element "i" is a NaN, then it
% will be preserved as NaN in the output.  At the edges of the matrix,
% where you cannot build a full rectangle, the values are unchanged. 
% 
% "matrixIn": original matrix
% "Nr": number of points used to smooth rows
% "Nc": number of points to smooth columns.  If not specified, Nc = Nr.
% 
% "matrixOut": smoothed version of original matrix
% 
% 
% 	Written by Greg Reeves, March 2009.
% 	Division of Biology
% 	Caltech
% 
% 	Inspired by "smooth2", written by Kelly Hilands, October 2004
% 	Applied Research Laboratory
% 	Penn State University
% 
% 	Developed from code written by Olof Liungman, 1997
% 	Dept. of Oceanography, Earth Sciences Centre
% 	G�teborg University, Sweden
% 	E-mail: olof.liungman@oce.gu.se

%
% Initial error statements and definitions
%
if nargin < 2, error('Not enough input arguments!'), end

N(1) = Nr; 
if nargin < 3, N(2) = N(1); else N(2) = Nc; end

if length(N(1)) ~= 1, error('Nr must be a scalar!'), end
if length(N(2)) ~= 1, error('Nc must be a scalar!'), end

%
% Building matrices that will compute running sums.  The left-matrix, eL,
% smooths along the rows.  The right-matrix, eR, smooths along the
% columns.  You end up replacing element "i" by the mean of a (2*Nr+1)-by- 
% (2*Nc+1) rectangle centered on element "i".
%

A=smooth2a(matrixIn,Nr,Nc);
matrixIn(1+Nr:end-Nr,1+Nc:end-Nc)=A(1+Nr:end-Nr,1+Nc:end-Nc);
matrixOut=matrixIn;




