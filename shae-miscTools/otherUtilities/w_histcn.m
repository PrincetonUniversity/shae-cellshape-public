function [count edges mid loc]=w_histcn(X,W,varargin)
%does the histcn function, but weighted

switch length(varargin)
    case 1
[count1 edges mid loc] = histcn(X, varargin{1});
    case 2
        [count1 edges mid loc] = histcn(X, varargin{1}, varargin{2});

    case 3
        [count1 edges mid loc] = histcn(X, varargin{1}, varargin{2},...
            varargin{3});

    case 4
        [count1 edges mid loc] = histcn(X, varargin{1}, varargin{2},...
            varargin{3}, varargin{4});
end





%siz= size(count1);
count=zeros(size(count1));
   W=W(all(loc,2))';
loc=loc(all(loc,2),:);
count=accumarray(loc,W);
%   locsub=sub3ind(siz,loc);
%count(sub3ind(siz,loc))=count(sub3ind(siz,loc))+W;

% for s=1:length(locsub)
% 
% count(locsub(s))=count(locsub(s))+W(s);
% %    
% %        if all(loc(s,:)>0)
% %           count(sub3ind(siz,loc(s,:)))=count(sub3ind(siz,loc(s,:)))+W(s);
% %        end
%        
% end
if size(count,2)~=1
for s=1:size(count,1)
    count(s,:)=count(s,:)/sum(count(s,:));
    
end
else
    count=count/sum(count);
end

if ~all(size(count)==size(count1))
    count=padarray(count,size(count1)-size(count),0,'post');
end
