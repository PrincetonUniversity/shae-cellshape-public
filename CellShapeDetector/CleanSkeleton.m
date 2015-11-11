function S=CleanSkeleton(bwin)


cc=bwconncomp(bwin);
maxArea=0;
for j=1:cc.NumObjects
    maxArea=max(maxArea,length(cc.PixelIdxList{j}));
end
for j=1:cc.NumObjects
    if length(cc.PixelIdxList{j})<maxArea
        bwin(cc.PixelIdxList{j})=0;
    end
    
end
bwin(1,:)=0;
bwin(:,1)=0;
bwin(end,:)=0;
bwin(:,end)=0;
i=1;
bwin_history(:,:,i)=bwin;
while sum(sum(bwmorph(bwin,'branchpoints')))~=0
    i=i+1;  
    bwin=bwmorph(bwin,'spur');
    bwin=bwmorph(bwin,'thin',Inf);

    bwin_history(:,:,i)=bwin;
%     imagesc(bwin+bwmorph(bwin,'branchpoints'));
 %pause(.1)

end

S=bwin_history(:,:,i);
for k=i-1:-1:1;
    ends=bwmorph(S,'endpoints');
    end_no=sum(ends(:));
    ends=bwmorph(ends,'thicken',1);
    ends=bwmorph(ends,'thin',1);
    pad=and(ends,((bwin_history(:,:,k)-bwin_history(:,:,k+1))>0));
    S=S+pad;
    
while end_no~=2
    K=find(bwmorph(S,'endpoints'));
    S(K(round(rand(1).*end_no+.5)))=0;
        ends=bwmorph(S,'endpoints');
    end_no=sum(ends(:));
end
    
    
 %   imagesc(S)
%pause(.1)
    
end
    ends=bwmorph(S,'endpoints');
    end_no=sum(ends(:));
    %if there are more than 2 end points)
while end_no~=2 && ~(end_no==1 &&  sum(S(:))==1);
    K=find(bwmorph(S,'endpoints'));
    if isempty(K)
        break
    end
    
    S(K(round(rand(1).*end_no+.5)))=0;
        ends=bwmorph(S,'endpoints');
    end_no=sum(ends(:));
    
end


S=bwmorph(S,'thin',Inf);
