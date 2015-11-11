function Imout=AreaFilter(Imin,minArea,maxArea,conn)
%takes a binary image Imin, and removes objects that are smaller than
%minArea and larger than maxArea with connectivity conn.

if nargin<4
    if ismatrix(Imin)==2
    conn=4;
    else
    conn=8;
    end
end
if nargin<3 
maxArea=Inf;
end

if isempty(maxArea);
    maxArea=Inf;
end

cc=bwconncomp(Imin,conn);
blobSizes=cellfun(@(x) length(x), cc.PixelIdxList);
cc.PixelIdxList(blobSizes<minArea |blobSizes>maxArea)=[];

cc.NumObjects=length(cc.PixelIdxList);

Imout=false(size(Imin));
Imout(cell2mat(cc.PixelIdxList'))=true;



