function varargout=weightedSmooth(varargin)
% Smooths 2D array data.
%
% function [mat1Out, mat2Out, ...] = smooth2a(xf, yf, zf, mat1In, mat2In, ... ,Nr,Nc)
% 
% This function smooths the data in matrixIn after weighting by the surface
% area of each rectangular face. The smooting uses a mean filter over a
% rectangle of size (2*Nr+1)-by-(2*Nc+1). It also takes into consideration
% the fact that the image has periodic boundary conditions along dimension
% 1 and that those periodic boundary conditions and the singularity at the 
% pole lead to a periodic boundary condition along the long axis with a
% phase shift of pi.

PERIODIC_REP_NUMBER = 3;

% check to see how the function is called
if nargin<4
     error('not enough inputs')
elseif nargin>=4
    % first three arguments are the x, y and z coordinates of the surface
    xf=varargin{1};
    yf=varargin{2};
    zf=varargin{3};
end

% default smoothing parameters if the last and second to last entries have
% a size larger than 1
if numel(varargin{end})~=1 && numel(varargin{end-1})~=1
    m=4;
    n=4;
else
    % use the last two parameters as the size scale for the smoothing
    m=varargin{end};
    n=varargin{end-1};
    varargin(end-1:end)=[];
end

% calculate the size of the matrix
[row,col]=size(xf);
% calculate the surface area of each element
[~,a]=surfArea(xf,yf,zf);
% remove the final row
% BPB: does this removal require that the first and last lines of xf
% be the same? If so, there should be a check to make sure that is true
a(end,:)=[];
% replicate for periodic boundary conditions
A=repmat(a,PERIODIC_REP_NUMBER,1);
% Along the long axis, the peroidic boundary conditions come back with a pi
% phase shift.
% BPB: Do these hard coded 20s come an expectation of the size of the
% matrix to be 40? Also, there should be a fliplr in here to account for
% the opposite direction
A=[circshift(A,-20),A,circshift(A,20)];
% smooth this surface area matrix using the same size scale along both
% dimensions
A=smooth2a(A,m,m);

% Number of outputs must be >=minargs and <=maxargs.
nargoutchk(1, length(varargin)-3);

% pull out the arrays to be weighted/smoothed. For each of them, weight
% them by the surface area and then smooth using the size scales m and n
for i=4:length(varargin)
    % use the name gc as the temporary holder of the data to be smoothed
    gc=varargin{i};
    % remove the last row if it has been replicated or is not the correct
    % size
    if size(gc,1)~=size(a,1)
        gc(end,:)=[];
    end
    
    % replicate for periodic boundary conditions
    Agc=repmat(a.*gc,PERIODIC_REP_NUMBER,1);
    % Along the long axis, the peroidic boundary conditions come back with a pi
    % phase shift.
    % BPB: Do these hard coded 20s come an expectation of the size of the
    % matrix to be 40?
    Agc=[circshift(Agc,-20),Agc,circshift(Agc,20)];
    % smooth by sizes m, n
    Agc=smooth2a(Agc,m,n);
    
    % remove the weighting -> back into a "density" metric
    Agc=Agc./A;
    
    % store the results
    varargout{i-3}=Agc(row:row+row-1,col+1:col+col);
    
end

