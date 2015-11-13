% Original copyright (c) 2009, Bruno Luong
% All rights reserved.
% Modifications by Benjamin Bratton, 2012, see notes below.
function A = convnfft(A, B, shape, dims, options)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
% 
%   C = CONVNFFT(A, B, SHAPE) controls the size of the answer C:
%       'full'   - (default) returns the full N-D convolution
%       'same'   - returns the central part of the convolution that
%                  is the same size as A.
%       'valid'  - returns only the part of the result that can be
%                  computed without assuming zero-padded arrays.
%                  size(C,k) = max([nak-max(0,nbk-1)],0).
%
%   C = CONVNFFT(..., SHAPE, DIMS) with DIMS is vector of dimensions where
%       the convolution will be carried out. By default DIMS is
%       [1:max(ndims(A),ndims(B))] (all dimensions). A and B must have the
%       same lengths on other dimensions.
%   C = CONVNFFT(..., SHAPE, DIMS, GPU)
%       GPU is boolean flag, see next
%
%   C = CONVNFFT(..., SHAPE, DIMS, OPTIONS)
%
%   OPTIONS is structure with following optional fields
%       - 'GPU', boolean. If GPU is TRUE Jacket/GPU FFT engine will be used
%       By default GPU is FALSE.
%       - 'Power2Flag', boolean. If it is TRUE, use FFT with length rounded
%       to the next power-two. It is faster but requires more memory.
%       Default value is TRUE.
%       - 'KernelFourier', boolean. If it is TRUE, the kernel is already in
%       the fourier domain.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
%   Author: Bruno Luong <brunoluong@yahoo.com>
%   History:
%       Original: 21-Jun-2009
%       23-Jun-2009: correct bug when ndims(A)<ndims(B)
%       02-Sep-2009: GPU/JACKET option
%       04-Sep-2009: options structure
%       16-Sep-2009: inplace product
%       16-Aug-2012: BPB: without calls to inplaceprod for use without calling C code
%       17-Aug-2012: BPB: options parameters includes k space kernel

if nargin<3 || isempty(shape)
    shape = 'full';
end

if nargin<5 || isempty(options)
    options = struct();
elseif ~isstruct(options) % GPU options
    options = struct('GPU', options);
end

nd = max(ndims(A),ndims(B));
% work on all dimensions by default
if nargin<4 || isempty(dims)
    dims = 1:nd;
end
dims = reshape(dims, 1, []); % row (needed for for-loop index)

% GPU enable flag
GPU = getoption(options, 'GPU', false);
% Check if Jacket is installed
GPU = GPU && ~isempty(which('ginfo'));

% fourier transform of the kernel is saved and reused
kSpaceKernel = getoption(options,'FourierKernel',false);

% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
switch lower(shape)
    % full size
    case 'full',
        ifun = @(m,n) 1:m+n-1;
    % same as input
    case 'same',
        ifun = @(m,n) ceil((n-1)/2)+(1:m);
    % only valid
    case 'valid',
        ifun = @(m,n) n:m;
    otherwise
        error('convnfft: unknown shape %s', shape);
end


classA = class(A);
classB = class(B);


if not(kSpaceKernel);
    ABreal = isreal(A) && isreal(B);
else
    ABreal = isreal(A);
    % warn if real-ness of kernel is ignored
    if and(ABreal,isreal(B))
        warning('convnfft:realnessOfKernel',...
            ['Treating inputs as real even though the k space kernel is real',...
            ' and therefore the real space kernel complex'])
    end
end

% Special case, empty convolution, try to follow MATLAB CONVN convention
if any(size(A)==0) || any(size(B)==0)
    szA = zeros(1,nd); szA(1:ndims(A))=size(A);
    szB = zeros(1,nd); szB(1:ndims(B))=size(B);
    % Matlab wants these:
    szA = max(szA,1); szB = max(szB,1);
    szC = szA;
    for dim=dims
        szC(dim) = length(ifun(szA(dim),szB(dim)));
    end
    A = zeros(szC,classA); % empty -> return zeros
    return
end

power2flag = getoption(options, 'Power2Flag', true);
if power2flag
    % faster FFT if the dimension is power of 2
    lfftfun = @(l) 2^nextpow2(l);
else
    % slower, but smaller temporary arrays
    lfftfun = @(l) l;
end

if GPU % GPU/Jacket FFT
    if kSpaceKernel
        error('convnfft:invalidCombination','k space kernels and GPU support cannot be used at the same time');
    end
    
    if strcmp(classA,'single')
        A = gsingle(A);
    else
        A = gdouble(A);
    end
    if strcmp(classB,'single')
        B = gsingle(B);
    else
        B = gdouble(B);
    end
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        % compute the FFT length
        if kSpaceKernel % if the kernel is in kspace, use its size
            l = size(B,dim);
        else
            m = size(A,dim);
            n = size(B,dim);
            l = lfftfun(m+n-1);
        end
        % We need to swap dimensions because GPU FFT works along the
        % first dimension
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
            B = permute(B, swap);
        end
        A = fft(A,l);
        B = fft(B,l);
        subs{dim} = ifun(m,n);
    end
else % Matlab FFT
    % Do the FFT
    subs(1:ndims(A)) = {':'};
    for dim=dims
        % compute the FFT length
        if kSpaceKernel
            m = size(A,dim);
            l = size(B,dim);
            n = l-m+1;
        else
            m = size(A,dim);
            n = size(B,dim);
            l = lfftfun(m+n-1);
            
        end
        % fft of input image
        A = fft(A,l,dim);
        % fft of convolution kernel
        if not(kSpaceKernel)
            B = fft(B,l,dim);
        end
        subs{dim} = ifun(m,n);
        
    end
end
 
if GPU
    A = A.*B;
    clear B
else
    % in place product using builtin matlab code
    A=A.*B;
end

% Back to the non-Fourier space
if GPU % GPU/Jacket FFT
    for dim=dims(end:-1:1) % reverse loop
        A = ifft(A,[]);
        % Swap back the dimensions
        if dim~=1 % do the work when only required
            swap = 1:nd;
            swap([1 dim]) = swap([dim 1]);
            A = permute(A, swap);
        end        
    end   
else % Matlab IFFT  
    A = ifftn(A);
end

% Truncate the results
if ABreal
    % Make sure the result is real
    A = real(A(subs{:}));
else
    A = A(subs{:});
end

% GPU/Jacket
if GPU
    % Cast the type back
    if strcmp(class(A),'gsingle')
        A = single(A);
    else
        A = double(A);
    end
end

end % convnfft


%% Get defaut option
function value = getoption(options, name, defaultvalue)
% function value = getoption(options, name, defaultvalue)
    value = defaultvalue;
    fields = fieldnames(options);
    found = strcmpi(name,fields);
    if any(found)
        i = find(found,1,'first');
        if ~isempty(options.(fields{i}))
            value = options.(fields{i});
        end
    end
end
