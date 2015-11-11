function varargout = trimDataEnds(varargin)


cut=varargin{1};
for i=1:nargin-1
    K=varargin{i+1};
   varargout{i}=K(1:end,cut+1:end-cut,:);
end
