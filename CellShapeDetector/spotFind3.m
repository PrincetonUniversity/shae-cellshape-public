function fe=spotFind3(Fin)
%finds spots on the unwrap cell contour. The output is the unwrap
%replicated 3x. The segmentation uses first performs a bandpass filter and
%then finds spots using 2 thresholds, 1 global and one moving. Circular
%boundary conditions are accounted for for spots that extend past the edge.
%end
gkernel=convnfft(gausswin(5),gausswin(5)');
globalthresh=.1;
Fin=normalizeRange(Fin);

%replicate 3x for boundary
Finm=[Fin;Fin;Fin];

Finm=convnfft(Finm,gkernel,'same');
Finm=imtophat(Finm,strel('disk',10));


%Finm=wiener2(Finm,3);
%repeat middle section 3x and use moving threshold
Finm=Finm(size(Fin,1)+1:2*size(Fin,1),:);
Finm=normalizeRange(Finm);


H=hessianMatrix([Finm;Finm;Finm],1);
Heig=hessianEig(H);
Htrace=sum(Heig,3);

FM=Htrace<0 & [Finm;Finm;Finm]>globalthresh;
%select only middle section
FM=repmat(FM(size(Finm,1)+1:2*size(Finm,1),:),[3,1]);

%pad ends, open, and remove sides, leaving top and bototm
FM=padarray(FM,[3,1],0,'both');
%FM=bwmorph(FM,'open');
FM=FM(4:end-3,:);
FM=FM(:,2:end-1);

%subIm=FM.*[Fin;Fin;Fin];
subIm=bwdist(~FM);
%subIm=windowNormalize(subIm,15);
    subIm=imhmax(subIm,.7);
    
 fshed=(watershed(-subIm)>0);
 FM(~fshed)=0;

%for structures that are completely circumferential, add back later
FZ2=bwlabel(FM);
FM2=(FZ2~=0 & FZ2==circshift(FZ2,[size(FZ2,1)/3,0]));

FM=bwmorph(FM,'fill');
FM=repmat(FM(size(Finm,1)+1:2*size(Finm,1),:),[3,1]);
FM=padarray(FM,[0,1],0,'both');
FM=imclearborder(FM,8);
FM=FM(:,2:end-1);

FM=FM+FM2;

FZ=bwlabel(FM,4);
%for very long structures that would get removed by imclearborder, add them
%back later
FZ1=FZ(1:size(Fin,1),:);
FZ2=FZ(2*size(Fin,1)+1:end,:);
FZ3=FZ1.*(FZ1==FZ2);

%using circshift, average the mask labels over whole rotations. sqrt makes 
% it so labels stay unique.
FZ=sqrt(FZ)*100;
FZ=round(FZ);
FZ=(FZ+circshift(FZ,[size(FZ,1)/3,0])+circshift(FZ,[-size(FZ,1)/3,0]));

FZ=padarray(FZ,[0,1],0,'both');

FZ=imclearborder(FZ,4);
FZ=FZ(:,2:end-1);

FZ=FZ+[FZ3;FZ3;FZ3];
regionLabels=unique(FZ(:));
regionLabels(regionLabels==0)=[];



Fin=[Fin;Fin;Fin];

if isempty(regionLabels)
    
    FZ=bwlabel(FM);
FZ=padarray(FZ,[0,1],0,'both');
FZ=(FZ+circshift(FZ,[size(FZ,1)/3,0])+circshift(FZ,[-size(FZ,1)/3,0]));
FZ=FZ(:,2:end-1);
FZ=FZ(1:size(Fin,1),:);
regionLabels=unique(FZ(:));
regionLabels(regionLabels==0)=[];
end

if isempty(regionLabels)
    G2=zeros(size(Fin,1),2*size(Fin,2));
    F2=zeros(size(Fin,1),2*size(Fin,2));
end


%%check subimages and split them if needed. 
for areas=1:length(regionLabels)
    
    
    reg = (FZ==regionLabels(areas));
    reg =bwlabel(reg);
    reg=(reg==min(max(reg(:)),2));
    
    FZ(FZ==regionLabels(areas) &~reg)=0;
end



%Area Filter

for areas=1:length(regionLabels);
    A=sum(sum(FZ==regionLabels(areas)));
    if A<20
        FZ(FZ==regionLabels(areas))=0;
        regionLabels(areas)=0;
    end
end

regionLabels(regionLabels==0)=[];
%FZ=imclearborder(FZ,8);
FZ=imfill(FZ,'holes');
fe=bwlabel(FZ,4);












