%%This is the main operating code for cell shape detector. It fits the
%3d stack of a cell with the contour that best describes the cell surface
%after convolution with PSF. 



%select image stack to load. Stacks must be ordered by Z then channel.
%if running "multirun", comment the first section.
%load the PSF
%largePSF

clc;
close all
clear all

%image folder
 
p1='/Users/jeffnguyen/Documents/Data/ExampleFittingProcess/'; 

p1='/Volumes/processedData/JeffreyNguyen/Simulated Images/FitTesting/';
cd(p1)
% stack_name ='TIMEPOINT001_20130824_162544-13';    %image name
% stack_name='TIMEPOINT001_20130703_225504-3';
stack_name='wavy20131016T030436';



ptail='conv3d3';    %string to append on end of output
p=strcat(p1,strrep(stack_name,'.tif',''));


%% Obtaining size of im


[row col] = size(imread([p,'.tif'], 'tif', 1));

image_para.nm_per_pixel=80; %AFM 62, Block 80old, Blocknew 75
image_para.Z_scale=.6;         %.6   %scaling between distance in the slide to distance moved in the image
image_para.stack_z_size=54;
image_para.stack_t_size=1;
image_para.Fstack_z_size=40;
image_para.Fstack_t_size=1;
image_para.stack_seperation_nm=80/.6;
image_para.Fstack_seperation_nm=100;

image_para.ZF_factor=...
    image_para.Z_scale*image_para.Fstack_seperation_nm/image_para.nm_per_pixel;
image_para.nm_per_zstack=...
image_para.Z_scale*image_para.stack_seperation_nm;

image_para.ZX_factor=image_para.nm_per_zstack/image_para.nm_per_pixel;

t_start=1;          %start and end time
t_end=1;
if t_end>image_para.stack_t_size
    t_end=image_para.stack_t_size;
end


pixel_interp=1;

%Flags
flag.fftconv=1;     %convolve with FFT (only 1)
flag.F1=0;                %is there an additional fluorescent image [0,1];
Fflag2=0;                  %is there a second fluorescent protein (only 0 for now)
flag.showfits=1;        %show 3d fits and protein fits
flag.gradient=0;         %use gradient vs max 

%Are poles in the image (only 1 for now)
flag.endpole=1;
flag.startpole=1;

flag.zplanar=1;         %is the cell roughly planar in Z
flag.straight=1;           %is the cell centerline roughly straight

flag.reverse=0;    %set 0 for shape then protein, 1 for protein then shape
%intensity_correction = -32768;
flag.psf4d=0;   %use full dvpsf? (only 0 for now)
intensity_correction = 0;
window=40;  %sizlosee of window size when reslicing cell in pixels
end_angles=20; %number of slices used to fit poles



if image_para.Fstack_t_size>1 && flag.F1==1
    mkdir(p)
end
%% Initialize cline paramaters

cline_para.iterations=700;
cline_para.group_iterations=40;%40
cline_para.it_loop=80;%70
cline_para.pts=40; %pts in circle (needs to be divisible by 4)
cline_para.iterations3d=300;
if rem(cline_para.pts,4)~=0;
    cline_para.pts=cline_para.pts+4-rem(cline_para.pts,4);
end
cline_para.radius=4*pixel_interp; %initial radius
cline_para.center=[pixel_interp*window+2,pixel_interp*window+2];
cline_para.threshold=.1;
cline_para.stiff=0.05; %.0005
cline_para.stiff3d=1; %.1
cline_para.alpha =.005 ;
cline_para.beta = .5;
% cline_para.alpha3d =00;
% cline_para.beta3d =500;
% cline_para.zalpha3d=0;
% cline_para.zbeta3d=500;
cline_para.alpha3d =00;
cline_para.beta3d =10;
cline_para.zalpha3d=0;
cline_para.zbeta3d=5;


cline_para.inter_frame_viscos_factor=0; %0
cline_para.kappa3d = .5; %1  %energy term
cline_para.kappa=10;
cline_para.gamma=10/pixel_interp; %step size
cline_para.gradflag=0;
cline_para.gradient_force =.5;%5
cline_para.interpolate_gradient_flag =0;
cline_para.show_flag=2; %1 to show end of every frame, 2 to show every compelted slice
cline_para.movie_flag=0;
cline_para.movie_title=[stack_name 'Allfit.avi'];



%% find and load metadata
%%% ----- filter shift patch

% find the last dash in the filename which separates images into individual cells
% this dash has ascii code 45
try
    dashPos = find(double(stack_name)==45,1,'last');
    metaFileName = stack_name;
    metaFileName(dashPos:end) = [];
    metaFileName = [metaFileName,'-Metadata.mat'];
    load(metaFileName);
    isMetaData = 1;
catch ME
    isMetaData = 0;
end

%%% ----- end filter shift patch

%% load raw hyper stack
%hyper_stackx = zeros(row, col);
hyper_stack=zeros(row,col,image_para.stack_z_size);
for k = 1:image_para.stack_z_size
    % determine which plane should be read in
    if flag.reverse==1 && flag.F1==1
        planeIndex = image_para.Fstack_t_size*image_para.Fstack_z_size+k;
    else
        planeIndex = k;
    end
    % read in the plane and store it
    hyper_stack(:,:,k) = double(imread([p,'.tif'], 'tif',...
            planeIndex)...
            +intensity_correction);
    
%  hyper_stack(:,:,k)=wiener2(hyper_stack(:,:,k));
%   if ~isempty(find(hyper_stack(:,:,k)<0, 1)) && k==1
%     H1=hyper_stack(:,:,k);
%     Hmask=(H1==0);
%     Hmask(smooth2a(Hmask,3,3)>0)=1;
%     
%   elseif ~isempty(find(hyper_stack(:,:,k)<0, 1))
%      H1=hyper_stack(:,:,k);
%     Hmin=min(H1(~Hmask));
%     H1(Hmask)=Hmin;
%     hyper_stack(:,:,k)=H1;
%   end
  

end
%%% ----- filter shift patch
if isMetaData
    % background subtract, shiftImgFilter, background add
    backgndTmp = quantile(hyper_stack(:),0.005);
    hyper_stack = hyper_stack-backgndTmp;
    for k2 = 1:image_para.stack_z_size
        if flag.reverse==1 && flag.F1==1
            planeIndex = image_para.Fstack_t_size*image_para.Fstack_z_size+k2;
        else
            planeIndex = k2;
        end
        
        filterString = totalImageInfo{planeIndex}.filterID;
        filterIdx = str2double(filterString(end-1));
        if isnan(filterIdx)
            filterIdx = 3; % no shifting of image
        end
    
        hyper_stack(:,:,k2) = shiftImgFilter(hyper_stack(:,:,k2),filterIdx);
        
    end
    hyper_stack = hyper_stack+backgndTmp;
end
%%% ----- end filter shift patch


hyper_stack=image_resize(hyper_stack,size(hyper_stack,1),...
    size(hyper_stack,2),round(image_para.ZX_factor*size(hyper_stack,3)));
image_para.imsize=size(hyper_stack);
hyper_stack=normalizeRange(hyper_stack);




%% Paramaterize centerline by distance instead of coordinate
z_level=.15;
[xyzs,zbw,zsq]=cellCenterline3d(hyper_stack,z_level,window,...
    image_para.imsize(3),flag);

%% Calculate tangent vectors of points and create new basis for detecting membrane
if size(xyzs,1)>2
    [Tv,Bv,Nv]=tbnVector(xyzs);
    %creates curvilinear coord around the centerline using TBN unit vectors,
    %B and N are rotated to be as parallel as possible to previous slice
else
    Tv=[1,0,0];
    Bv=[0,1,0];
    Nv=[0,0,1];
end


plane_num=size(Tv,1);
%make the first and last 'endround' tbn vectors the same so nothing
%strange happens at the ends.
if plane_num>10
    endround=5;
else
    endround=round(plane_num/2);
end

for i=1:endround
    Tv(plane_num-i+1,:)=Tv(plane_num-endround+1,:);
    Bv(plane_num-i+1,:)=Bv(plane_num-endround+1,:);
    Nv(plane_num-i+1,:)=Nv(plane_num-endround+1,:);
    Tv(i,:)=Tv(endround,:);
    Bv(i,:)=Bv(endround,:);
    Nv(i,:)=Nv(endround,:);
end
%create a 2*window +1 square around each point for interpolationg using
%the B and N vectors

xslice=zeros(pixel_interp*2*window+1,...
    pixel_interp*2*window+1,size(Bv,1));
yslice=xslice;
zslice=xslice;

[J,K]=meshgrid(-window:1/pixel_interp:window,...
    -window:1/pixel_interp:window);
        for i=1:size(Bv,1)
            
            
            xslice(:,:,i)=J(:,:)*Nv(i,1)+K(:,:)*Bv(i,1)+xyzs(i,1);
            yslice(:,:,i)=J(:,:)*Nv(i,2)+K(:,:)*Bv(i,2)+xyzs(i,2);
            zslice(:,:,i)=J(:,:)*Nv(i,3)+K(:,:)*Bv(i,3)+xyzs(i,3);
            
        end
        
        



%use points to interpolate, XY in matlab is messed me up.. but this
%works
%if using gradient, convolve stack with sobel operator in each
%direction and find magnitude

if cline_para.gradflag==0
    
    V=interp3(hyper_stack,yslice,xslice,zslice,'*linear',0);
    
elseif cline_para.gradflag==1
    kernal=fspecial('sobel');
    kernal3(:,:,1)=kernal;
    kernal3(:,:,2)=2*kernal;
    kernal3(:,:,3)=kernal;
    kernal3r=kernal3;
    kernal3p=permute(kernal3,[2,3,1]);
    kernal3c=permute(kernal3,[3,1,2]);
    hostack=convn(hyper_stack,kernal3r,'same');
    vertstack=convn(hyper_stack,kernal3c,'same');
    pagestack=convn(hyper_stack,kernal3p,'same');
pagestack(:,:,1)=0;
pagestack(:,:,end)=0;

    grad_stack=sqrt(hostack.^2+vertstack.^2+pagestack.^2);
    V=interp3(grad_stack,yslice,xslice,zslice,'linear',0);
    %V=smooth3(V);
    
    clear hostack vertstack pagestack
end
if flag.zplanar==1 && flag.straight==1
    %V=smooth3(V,'box',[1,1,3]);
    %        V=smooth3(V,'gaussian',[5,5,5],3*pixel_interp);
end

if flag.showfits~=1
    clear xlice yslice zslice
end



%% Extract boundaries cell in each slice and centers at the end for poles
%this works in the xyzslice coordinate system

[xys,c1,c2] = f_CircDet8(V, cline_para); %gets center line and end pts
close all
c1=c1-[window*pixel_interp+1,window*pixel_interp+1];
c2=c2-[window*pixel_interp+1,window*pixel_interp+1];

firstloop=xys(:,:,1)-repmat(c1,length(xys(:,:,1)),1);
lastloop=xys(:,:,plane_num)-repmat(c2,length(xys(:,:,1)),1);


c1=c1/pixel_interp;
c2=c2/pixel_interp;


%% Transform back to image space, smooth'
%substract center offset from fitting in window
xys_all = xys-repmat(pixel_interp*window+1,size(xys));

xys_final=zeros([size(xys_all,1),size(xys_all,2)+1,size(xys_all,3)]);

for i=1:size(xys_all,3)
    %mutliply each slice by basis vectors to return to real space
    xyz= xys_all(:,:,i)*[Nv(i,:)',Bv(i,:)']'/pixel_interp;
    xys_final(:,:,i)=repmat(xyzs(i,:),[size(xys_all,1),1,1])+xyz;
end

xys_final(:,1:2,:)=xys_final(:,1:2,:);
xys_final(:,3,:)=xys_final(:,3,:);
if flag.showfits~=1
    clear V
end
%% POlES
if flag.gradient==0
[P_start,P_end,theta_start,theta_end, Bp_start,Bp_end,start_pt,end_pt]=...
    cellPoleCoordSys(hyper_stack,Nv,Bv,c1,c2,xyzs,pixel_interp,window,flag,end_angles);    
elseif flag.gradient==1
[P_start,P_end,theta_start,theta_end, Bp_start,Bp_end,start_pt,end_pt]=...
    cellPoleCoordSys(grad_stack,Nv,Bv,c1,c2,xyzs,pixel_interp,window,flag,end_angles);
clear grad_stack
end

%% Extract Pole Boundaries from P stack
%cline_para.stiff=.001;

cline_para.group_iterations=00;
%cline_para.iterations=500;


if flag.startpole==1
    
    xys_end_s = f_CircDet3(P_start, cline_para,firstloop);
    xys_end1=xys_end_s-repmat(pixel_interp*window+1,size(xys_end_s));
end
if flag.endpole==1
    xys_end_e = f_CircDet3(P_end, cline_para,lastloop);
    xys_end2=xys_end_e-repmat(pixel_interp*window+1,size(xys_end_e));
    
end


close all


%% Transform Poles back to image space
clear P_start P_end
[xpe,ype,zpe,l_end]=...
    cellPoleCoordSys2Image(xys_end2,Nv(end,:),Bp_end,cline_para,end_angles,pixel_interp);
[xps,yps,zps,l_start]=...
    cellPoleCoordSys2Image(xys_end1,Nv(1,:),Bp_start,cline_para,end_angles,pixel_interp);



%% Fix Ends and reparameterize into coordinate system that matches the body
%This converts pole points into spherical coordinates, rotates the
%coordinate axis to a new phi and theta, and then fixes the angle to
%replace it at the poles. It also calculates curvatures in the original
%coordinate system to use later for replacing the singularity at the
%the tips of the cell.
[xp2s,yp2s,zp2s,start_slices]=PoleRotate(l_start,cline_para,-Tv(1,:),...
    xps,yps,zps,theta_start);
xp2s=xp2s+start_pt(1);
yp2s=yp2s+start_pt(2);
zp2s=zp2s+start_pt(3);
[xp2e,yp2e,zp2e,end_slices]=PoleRotate(l_end,cline_para,Tv(end,:),...
    xpe,ype,zpe,theta_end);
xp2e=xp2e+end_pt(1);
yp2e=yp2e+end_pt(2);
zp2e=zp2e+end_pt(3);


%% Combine poles with body
skip=1;  %skip every few points to clean up


x0=squeeze(xys_final(:,1,1:skip:end));
y0=squeeze(xys_final(:,2,1:skip:end));
z0=squeeze(xys_final(:,3,1:skip:end));
if size(xys,3)>=2
    xys(:,:,1)=[];  %first and last rings are identical to poles
    xys(:,:,end)=[];
elseif size(xys,3)==1
    xys(:,:,1)=[];
end

%cut off ends and attach poles in proper orientation.
xp2s=flipdim(xp2s,1);
yp2s=flipdim(yp2s,1);
zp2s=flipdim(zp2s,1);
% 
xp2e=flipdim(xp2e,2);
yp2e=flipdim(yp2e,2);
zp2e=flipdim(zp2e,2);
xp2s=flipdim(xp2s,2);
yp2s=flipdim(yp2s,2);
zp2s=flipdim(zp2s,2);

xp2e=circshift(xp2e,[0,-1]);
yp2e=circshift(yp2e,[0,-1]);
zp2e=circshift(zp2e,[0,-1]);

xp2s=circshift(xp2s,[0,-1]);
yp2s=circshift(yp2s,[0,-1]);
zp2s=circshift(zp2s,[0,-1]);
x0=x0(:,2:end);
y0=y0(:,2:end);
z0=z0(:,2:end);



if ~isempty(x0)
x1=[xp2s',x0,xp2e'];
y1=[yp2s',y0,yp2e'];
z1=[zp2s',z0,zp2e'];
else
%     xp2e=xp2e(2:end,:);
%     yp2e=yp2e(2:end,:);
%     zp2e=zp2e(2:end,:);

    x1=[xp2s',x0,xp2e'];
y1=[yp2s',y0,yp2e'];
z1=[zp2s',z0,zp2e'];
end

        
%% estimate centerline
if size(x0,2)>0
    centerline=zeros(size(x0,2),3);
    centerline(:,1)=mean(x0);
    centerline(:,2)=mean(y0);
    centerline(:,3)=mean(z0);
else
    centerline=[];
end

cl=zeros(size(x1,2),3);
cl(:,1)=mean(x1);
cl(:,2)=mean(y1);
cl(:,3)=mean(z1);



%% Load PSF
upFactor=2;
CLscale=4;
[x,y,z,eg]=activeSurfaceFit(x1,y1,z1,hyper_stack,cline_para,image_para,upFactor,flag);

        membrane2=coord2image3d(x,y,z,image_para.imsize,CLscale,1);
        mfill=imfill(double(membrane2~=0),'holes');
        CL=[mean(x)' mean(y)' mean(z)'];
        M3=bwdist(~mfill).^2;
        M3=M3/max(M3(:))*CLscale;
        CL = ActiveContourFit(M3, cline_para, CL*CLscale);
        
   %% Create Fluorescent map by interpolating

if flag.F1==1
    xf=x;
    yf=y;
    zf=z;
    [x0,y0,z0]=...
    trimDataEnds(end_slices,xf,yf,zf);

    
    F=zeros(row,col,image_para.Fstack_z_size);
    for tt=1:image_para.Fstack_t_size
        % if theres a Fluorescent stack, open it and scale it so that the
        %axis are equal and have its center be the same as the membrane stack
        F1=zeros(row,col,image_para.Fstack_z_size);
        for k = 1:image_para.Fstack_z_size
            if flag.reverse==0
                planeIndex =  image_para.stack_z_size*image_para.stack_t_size+(tt-1)*image_para.Fstack_z_size+k;
                
            else
                planeIndex = (tt-1)*image_para.Fstack_z_size+k;
                
            end
            F1(:,:,k) = double(imread([p,'.tif'], 'tif',planeIndex))+...
                    intensity_correction;
        end
        %%% ----- filter shift patch
        % check for filter information
        if isMetaData
            % background subtract, shiftImgFilter, background add
            backgndTmp = quantile(F1(:),0.005);
            F1 = F1-backgndTmp;
            for k2 = 1:image_para.Fstack_z_size
            if flag.reverse==0
                planeIndex =  image_para.stack_z_size*image_para.stack_t_size+(tt-1)*image_para.Fstack_z_size+k2;
                
            else
                planeIndex = (tt-1)*image_para.Fstack_z_size+k2;
                
            end

                filterString = totalImageInfo{planeIndex}.filterID;
                filterIdx = str2double(filterString(end-1));
                if isnan(filterIdx)
                    filterIdx = 3; % no shifting of image
                end

                F1(:,:,k2) = shiftImgFilter(F1(:,:,k2),filterIdx);

            end
            F1 = F1+backgndTmp;
        end 
        %%% ----- end filter shift patch
        F1=image_resize(F1,size(F1,1),...
    size(F1,2),round(image_para.ZF_factor*size(F1,3)));
        
        % This will show the fit outlines on the fluorescent channel
        if flag.showfits==1
            FV=interp3(F1,yslice,xslice,zslice,'linear*',0);
            for k=1:size(xys,3)
                subplot(1,2,1)
                imagesc(squeeze(FV(:,:,k)))
                hold on
                plot(xys(:,1,k), xys(:,2,k), 'g.-');
                
                axis equal
                subplot(1,2,2)
                imagesc(squeeze(V(:,:,k)))
                hold on
                plot(xys(:,1,k), xys(:,2,k), 'g.-');
                axis equal
                pause(.1)
            end
        end
        
        %This is set up to measure the intensities at the membrane for each
        %point.
        Ff=interp3(F1,yf, xf,zf,'*linear');
        Ff0=interp3(F1,y0, x0,z0,'*linear');
        I=interp3(hyper_stack,yf, xf,zf,'*linear');
        
        %max radial projection
        
        Ffr=zeros(size(x1));
        r_steps=0.5:.05:2;  %radius range
        rkx=zeros(length(r_steps),1);
        rky=rkx;
        rkz=rky;
        
        for i=1:size(x1,2)
            for k=1:size(x1,1)
                
                rkx=round(((xf(k,i)-CL(i,1))*r_steps+CL(i,1)));
                rky=round(((yf(k,i)-CL(i,2))*r_steps+CL(i,2)));
                rkz=round(((zf(k,i)-CL(i,3))*r_steps+CL(i,3)));
                
                rkx(rkx<1)=1;
                rky(rky<1)=1;
                rkz(rkz<1)=1;
                
                rkx(rkx>row)=row;
                rky(rky>col)=col;
                rkz(rkz>size(F1,3))=size(F1,3);
                
                for ii=1:length(r_steps)
                    F2(ii)=F1(rkx(ii),rky(ii),rkz(ii));
                    
                end
                Ffr(k,i)=max(F2);
                %Ftest(k,i,:)=F2;
                
            end
        end
        
        clear F1
    end
end
     
 %%       
        
x=x*image_para.nm_per_pixel;
y=y*image_para.nm_per_pixel;
z=z*image_para.nm_per_pixel;
CL=CL*image_para.nm_per_pixel/CLscale;


%% Calculate curvature and centerline
%gm and gc are mean and gaussian curvatures, xfyfzf are same as xyz,
%but closed at the end so that top and bottom values repeat(for
%visualization of closed cell.


start_pt2=start_pt*image_para.nm_per_pixel;
end_pt2=end_pt*image_para.nm_per_pixel;
[xf,yf,zf,gc,gm,kplus,kminus]=cellCurvature(x,y,z,Nv,Bv);
[r,CL]=CLRadiusFix(xf,yf,zf,CL);
[x0,y0,z0,gc0,gm0,kplus0,kminus0,r0]=...
    trimDataEnds(end_slices,xf,yf,zf,gc,gm,kplus,kminus,r);

if flag.F1==1
    if size(Ff)<size(xf,1)
        Ff=[Ff;Ff(1,:)];
        Ffr=[Ffr;Ffr(1,:)];
    end
    [Ffr0,Ff0,I0]=trimDataEnds(end_slices,Ff,Ffr,I);
end

%% Calculate surface area, volume and centerline curvature;
Area=surfArea(xf,yf,zf);
Volume=surfVolume(xf,yf,zf);
K=Curvature(CL,10);
T=Twist(CL,10);




%% Save files
if flag.F1==1
    if image_para.Fstack_t_size>1
        p2=[p '\' num2str(tt,'%02d')];
    else
        p2=[p,ptail];
    end
    
    save(p2,'p2','xf','yf','zf','x0','y0','z0','Ff','Ffr', 'Ff0','gc','gc0','gm0','kplus0','kminus0'...
        ,'r','r0','gm','kplus','kminus','CL','Area','Volume','K','T','cline_para','image_para','eg','I')
    
    
end

if flag.F1==0
    p2=[p,ptail];
    
    save(p2,'p2','xf','yf','zf','x0','y0','z0','gc','gc0','gm0','kplus0','kminus0'...
        ,'r','r0','gm','kplus','kminus','CL','Area','Volume','K','T','cline_para','image_para','eg')
    
end



%       clear xpe ype zpe xpef ypef zpef xp2ef yp2ef zp2ef xps yps zps xp2sf yp2sf zp2sf xpsf ypsf zpsf
%       clear xp2 yp2 zp2 xp2e yp2e zp2e xp2s yp2s zp2s Ftest cl centerline

%%

if Fflag2==1 %%this is for DNA stain or other signal on interior of cell
    warning('cellShape:untestedParameter:Fflag2','Fflag2 set to true. This is an untested mode.');
    F=zeros(row,col,image_para.Fstack_z_size);
    
    for k = 1:Fstack_z_size
        F(:,:,k) = double(imread(stack_name, 'tif',image_para.stack_z_size*image_para.stack_t_size*2+k)+intensity_correction);
    end
    F=image_resize(F,size(F,1),size(F,2),round(image_para.ZF_factor*size(F,3)));
    [X,Y,Z]=meshgrid(1:col,1:row,1:image_para.stack_z_size);
    [FX,FY,FZ]=meshgrid(1:col,1:row,1:size(F,3));
    F1=interp3(FX,FY,FZ-size(F,3)/2,F,X,Y,Z-image_para.stack_z_size/2,'nearest');
    
    for k=1:image_para.stack_z_size %repeat the ends of the fluorescent stack so all points beyond the end will look like the end. This is for extrapolation
        if k<image_para.stack_z_size/2 && all(all(isnan(F1(:,:,k))))
            F1(:,:,k)=F(:,:,1);
        end
        if k>image_para.stack_z_size/2 && all(all(isnan(F1(:,:,k))))
            F1(:,:,k)=F(:,:,end);
        end
    end
    clear F
    %This is set up to measure the average F intensity inside the
    %contour of each slice
    
    Ff2=zeros(size(xys,3),1);
    for i=1:size(xys,3)
        F=interp3(X,Y,Z,F1,yslice(:,:,i),xslice(:,:,i),zslice(:,:,i),'linear',NaN);
        bwmask=poly2mask(xys(:,1,i),xys(:,2,i),2*window+1,2*window+1);
        bwmask(isnan(F))=0;
        Ff2(i)=mean(F(bwmask));
        %             imagesc(F(:,:,i));axis equal
        %             hold on
        %             plot(xys(:,1,i),xys(:,2,i))
        %             pause(.1)
        %
    end
    
    
    %             Ff2=interp3(X,Y,Z,F1,centerline(:,2)/nm_per_pixel,...
    %                 centerline(:,1)/nm_per_pixel,centerline(:,3)/nm_per_pixel);
    clear F FX FY FZ
    
    
else Ff2=[];
end





%% Save the relevent variables


disp('done');

%% plot
%surf(x0(1:20,90:end),y0(1:20,90:end),z0(1:20,90:end));
%
% surf(xf,yf,zf,gc)
% axis equal
% caxis([min(min(gc))/5,-min(min(gc))/5]);
