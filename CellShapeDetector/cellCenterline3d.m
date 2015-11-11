function [xyzs,zbw,zsq,r_init]=cellCenterline3d(hyper_stack,z_level,window, stack_z_size2,flag)
%%
%            stack_z_size2=image_para.imsize(3);
%%
upFactor=2;
%threshold levels for Z and Y projections, adjust these to get good
%centerline
% hyper_stack2=smooth3(hyper_stack,'box',[3,3,3]);

se = strel('disk',100);
%z project, normalize, and filter image
% if flag.zplanar
% zstd=std(reshape(hyper_stack,[],stack_z_size2));
% zsq=hyper_stack(:,:,find(zstd==max(zstd),1,'first'));
% else
zsq=squeeze(mean(hyper_stack,3)); %use mean projection
%end
%%
zsq=interp2(zsq,upFactor,'*linear');       %interpolate image
zsq=imtophat(zsq,se);
zsq=normalizeRange(zsq);            %Normalize between 0 and 1
%%
%filter
zsq=smooth2a(zsq,3,3);
%zsq=medfilt2(zsq,[10,10]);
%%
zsq_bw=im2bw(zsq,z_level);
counter=1;

while sum(sum(imclearborder(zsq_bw)))<sum(zsq_bw(:))/2 && counter<8
    zsq_bw=im2bw(zsq,z_level+counter*.05);
    counter=counter+1;
end



zsq_bw=bwmorph(zsq_bw,'dilate');
zsq_bw=imfill(zsq_bw,'holes');
zsq_bw=bwmorph(zsq_bw,'erode');
zsq_bw([1,end],:)=0;
zsq_bw(:,[1,end])=0;
%%
zsq_bw=im2bw(smooth2a(double(bwdist(~zsq_bw)/max(max(bwdist(~zsq_bw))))...
    ,3,3),.5);

% skeletonize and removes spurious chains from 'thin' to get longest chain

zbw=bwmorph(zsq_bw,'thin',Inf);
zbw=CleanSkeleton(zbw);
zbw=bwmorph(zbw,'thin',Inf);
xyi=SkelOrder(zbw);  %convert BW image to ordered coordinates
%xyi2=unique(round(xyi/4-3/4),'rows');
%%
%xyi=xyi/2^upFactor-(2^upFactor-1)/2^upFactor; %divide and subtract to fix initial interpolation

%SvT changed the signs
xyi=(xyi-1)/2^upFactor+1; %divide and subtract to fix initial interpolation


%%

%take points in XY and find the Z coordinates cooresponding with
%centerline, ysq is the slice through the stack at the skeleton pts

se = strel('disk',5);
ysq=zeros(stack_z_size2,length(xyi));
for i = 1:size(xyi,1)
    ysq(:,i)=hyper_stack(round(xyi(i,1)),round(xyi(i,2)),:);
end

%ysq=imtophat(ysq,se);
ysq=(ysq-min(min(ysq)))/max(max(ysq-min(min(ysq))));

ysq=smooth2a(ysq,3,3);
ysq=wiener2(ysq,3);
ysq=normalizeRange(ysq);
ysq=bsxfun(@rdivide,ysq,max(ysq));    
ysq(ysq<.5)=0;
%imagesc(ysq)
ysq=wiener2(ysq,3);
yplot=zeros(1,length(xyi));
ymax = yplot;
yup=yplot;
ydown=yplot;
if flag.straight==1
    ysq=smooth2a(ysq,1,100);
end

%     estimate the cell center by finding the mean between the top half peak
%     and the bottom half peak (top and bottom of cell membranes)

ydown(1)=find(ysq(1:round(stack_z_size2/2),1)==max(ysq(1:round(stack_z_size2/2),1)),1,'first');
yup(1)=find(ysq(round(stack_z_size2/2):end,1)==max(ysq(round(stack_z_size2/2):end,1)),1,'first')+round(stack_z_size2/2);
yplot(1)=(yup(1)+ydown(1))/2;
for i = 2:length(xyi)
    lower=max(round(yplot(i-1))-window,1):(round(yplot(i-1)));
    upper=(round(yplot(i-1))):min(round(yplot(i-1))+window,stack_z_size2);
    if isempty(lower)
        lower=ydown(i-1);
    end
    if isempty(upper)
        upper=yup(i-1);
    end
    ydown(i)=mean(find(ysq(lower,i)==max(ysq(lower,i))));
    yup(i)=mean(find(ysq(upper,i)==max(ysq(upper,i))))+min(upper);
    
    yplot(i)=(yup(i)+ydown(i))/2;
end

% if isempty(lower) || isempty(upper)
for i=1:size(ysq,2)
    yplot(i)=sum(ysq(:,i)'.*(1:length(ysq(:,i))))/sum(ysq(:,i));
end
% end

yplot=smooth(yplot,30);

if flag.zplanar==1 % if cell is approximately planar, just find mean height
    yplot=ones(size(yplot))*mean(yplot);
end

if flag.straight==1 && size(xyi,1)~=1
    xyi(2,:)=xyi(end,:);
    xyi(3:end,:)=[];
end

% xyzi has points of centerline
xyzi=zeros(size(xyi,1),3);
for i=1:size(xyi,1)
    xyzi(i,3)=yplot(i);
    xyzi(i,1)=xyi(i,1);
    xyzi(i,2)=xyi(i,2);
end

% xyzi=flipdim(xyzi,1);
% trim head and tail if needed
%trim=5;
parse=1;
if length(xyzi)<parse
    parse=1;
end

%     xyzi(end-trim+1:end,:)=[];
%     xyzi(1:trim,:)=[];
xyzi=xyzi([1:parse:end-parse,end],:);
%  clear  zbw zsq ysq
if sum(ysq(:,1))<sum(ysq(:,end))
    xyzi=flipdim(xyzi,1);
end


if size(xyzi,1)~=1
    s=[0,cumsum(sqrt(diff(xyzi(:,1)).^2+diff(xyzi(:,2)).^2+diff(xyzi(:,3)).^2))']';
    xyzs=zeros(floor(max(s))+1,3);
    xyzs(:,1)=interp1(s,xyzi(:,1),0:1:max(s),'spline');
    xyzs(:,2)=interp1(s,xyzi(:,2),0:1:max(s),'spline');
    xyzs(:,3)=interp1(s,xyzi(:,3),0:1:max(s),'spline');
    
else
    xyzs=xyzi;
end
X=xyzs(:,1);Y=xyzs(:,2);Z=xyzs(:,3);
t=1:length(X);

X=fit(t',X,'smoothingSpline','SmoothingParam', .3);
X=X(t);
Y=fit(t',Y,'smoothingSpline','SmoothingParam', .3);
Y=Y(t);
Z=fit(t',Z,'smoothingSpline','SmoothingParam', .3);
Z=Z(t);

xyzs=[X,Y,Z];
trim=0;

xyzs(end-trim+1:end,:)=[];
xyzs(1:trim,:)=[];
xyzs=xyzs(1:end,:);
xyzs=smooth2b(xyzs,4,0);
zbw_d=bwdist(~zsq_bw);
r_init=double(mean(zbw_d(zbw))/(2^upFactor));