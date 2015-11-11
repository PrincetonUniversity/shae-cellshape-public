%%This is the main operating code for cell shape detector. It fits the
%3d stack of a cell with the contour that best describes the cell surface
%after convolution with PSF.
function Cell_shape_detectorPolymerFit(folderName,folderIdx,fileGroupSize);
% %%rm
% if nargin<1
% % image folder
%  p1='/Volumes/processedData/JeffreyNguyen/20130117_M63 FB83 Images/';
% %  cd(p1)
%  stack_name ='130117_183154timeseries-8';    %image name
% else
%     [p1,stack_name,~] = fileparts(fileName);
% end
% %%

% file group size is how many files to group together in one run of the
% main code
if nargin<3
    fileGroupSize = 1;
end

if fileGroupSize>1
    if numel(folderIdx)>1
        warning('cellShape:invalidParams','Suppressing extra folderIdx''s and using fileGroupSize');
    end
    folderIdx = folderIdx(1):folderIdx(1)+fileGroupSize-1;
    
end

disp(['Processing : ',folderName,'. Files : ',num2str(folderIdx),'.']);


tic;
[image_para,flag,cline_para] = paramsText(fullfile(folderName,'cell_shape_settings.txt'));
disp(['Processing params: ',num2str(toc),' seconds.']);
%%


% string to append on end of output
ptail='conv3d3';

fileList = dir([folderName,filesep,'*.tif']);


%% loop over files ------

for kkFile = 1:numel(folderIdx);
    tic;
    p1 = folderName;
    
    clearvars -except image_para flag cline_para ptail folderName folderIdx kkFile fileList p1
    stack_name = fileList(folderIdx(kkFile)).name;
    [~,stack_name,~] = fileparts(stack_name);
%    try
        
        % full file path to name of image
        p=strcat(p1,filesep,strrep(stack_name,'.tif',''));
        
        % ---- Parameters that are not used within structures
        pixel_interp = flag.pixel_interp;
        intensity_correction = flag.intensity_correction;
        window = flag.window;
        end_angles = flag.end_angles;
        pixel_interp = flag.pixel_interp;
        Fflag2 = flag.Fflag2;
        % ----
        
        % ---- Metadata for the image, '-' is 45
        p_metaIdx = find(double(p1),45,1,'last');
        p_meta = [p1(1:p_metaIdx-1),'metadata.mat'];
        try
            clear totalImageInfo
            load(p_meta);
            isMetadata = true;
        catch ME
            warning('CellShapeDetector:missingMetaData',...
                'Metadata.mat file missing.');
            isMetadata = false;
        end
            
        % ----
        
        %% Obtaining size of im
        [row col] = size(imread([p,'.tif'], 'tif', 1));
        
        % If there are multiple time points and there is a fluorescent channel to
        % be mapped, get a folder ready for them
        %if image_para.Fstack_t_size>1 && flag.F1==1
        %    mkdir(p)
        %end
        
        %% load raw hyper stack
        
        %hyper_stackx = zeros(row, col);
        
        if ~flag.F1
            image_para.stack_z_size=min(length(imfinfo([p,'.tif'])),...
                image_para.stack_z_size);
        end
        
        hyper_stack=zeros(row,col,image_para.stack_z_size);
        for k = 1:image_para.stack_z_size
            if flag.reverse==1 && flag.F1==1
                tempImageIdx = image_para.Fstack_t_size*image_para.Fstack_z_size+k;      
            else
                tempImageIdx = k;   
            end
            % read in the appropriate plane of the tif file
            tempImage = double(imread([p,'.tif'], 'tif',tempImageIdx)) + intensity_correction;
            
            % --- shift images if the filter is known
            if isMetadata
                filterUsed = totalImageInfo{tempImageIdx}.filterID;
                % current structure puts the filter number at
                % the end, such as 641/75 nm (3) is filter 3
                filterUsed = str2double(filterUsed(end-1));
                 hyper_stack(:,:,k) = shiftImgFilter(tempImage,filterUsed);
            else
                 hyper_stack(:,:,k) = tempImage;
            end
            % ---
            
            %  hyper_stack(:,:,k)=wiener2(hyper_stack(:,:,k));
            %             if ~isempty(find(hyper_stack(:,:,k)<0, 1)) && k==1
            %                 H1=hyper_stack(:,:,k);
            %                 Hmask=(H1==0);
            %                 Hmask(smooth2a(Hmask,3,3)>0)=1;
            %
            %             elseif ~isempty(find(hyper_stack(:,:,k)<0, 1))
            %                 H1=hyper_stack(:,:,k);
            %                 Hmask=(H1==0);
            %                 Hmin=min(H1(~Hmask));
            %                 H1(Hmask)=Hmin;
            %                 hyper_stack(:,:,k)=H1;
            %             end
            %
            
        end
        hyper_stack=image_resize(hyper_stack,size(hyper_stack,1),...
            size(hyper_stack,2),round(image_para.ZX_factor*size(hyper_stack,3)));
        image_para.imsize=size(hyper_stack);
        
        %pedistal subtraction
        pedMask=false(image_para.imsize);
        pedMask(1:3,:,:)=true;
        pedMask(:,1:3,:)=true;
        pedMask(end-2:end,:,:)=true;
        pedMask(:,end-2:end,:)=true;
       % pedMask(hyper_stack==0)=false;
        
        pedistal=median(hyper_stack(pedMask & (hyper_stack~=0)));
        hyper_stack(hyper_stack<pedistal)=pedistal;
        hyper_stack=normalizeRange(hyper_stack);
        if cline_para.gradflag
            hyper_stack=hyper_stack/max(max(max(smooth3(hyper_stack,'box',[3,3,3])))); %smooth abit so max is a smoothed max
        end
        
        
        %% Paramaterize centerline by distance instead of coordinate
        %         z_level=.2;
        
        
        % Ben adds for his centerline detection method, works pretty well
        % for starting off in triangular meshes as well
        if isfield(flag,'bpbFlag2')
        else
            flag.bpbFlag2  = false;
        end
        
        % use some centerline method
        if flag.bpbFlag2
            [xyzs,zbw,zsq]=cellCenterline3_bpb(hyper_stack,flag.z_level,window,...
                image_para.imsize(3),flag);
        else
            [xyzs,zbw,zsq]=cellCenterline3d(hyper_stack,flag.z_level,window,...
                image_para.imsize(3),flag);
        end
        
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
        
        if ~cline_para.gradflag
            
            V=interp3(hyper_stack,yslice,xslice,zslice,'*linear',0);
            
        else
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
            V(1,:,:)=0;
            V(end,:,:)=0;
            V(:,1,:)=0;
            V(:,end,:)=0;
            % BPB notes, first and last plane shouldn't be zero
            %V(:,:,1)=0;
            %V(:,:,end)=0;
            
            
            
            
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
        disp(['Finished first round of fits in ',stack_name,'. ',...
            'Elapsed time on this cell: ',num2str(toc)]);
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
        
        %% POlES
        if cline_para.gradflag==0
            [P_start,P_end,theta_start,theta_end, Bp_start,Bp_end,start_pt,end_pt]=...
                cellPoleCoordSys(hyper_stack,Nv,Bv,c1,c2,xyzs,pixel_interp,window,flag,end_angles);
        elseif cline_para.gradflag==1
            [P_start,P_end,theta_start,theta_end, Bp_start,Bp_end,start_pt,end_pt]=...
                cellPoleCoordSys(grad_stack,Nv,Bv,c1,c2,xyzs,pixel_interp,window,flag,end_angles);
        end
        
        
        %% Extract Pole Boundaries from P stack
        %cline_para.stiff=.001;
        
        cline_para.beta=cline_para.beta/10;
        
        
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
        CLscale=4;
        [x,y,z,eg]=activeSurfaceFit(x1,y1,z1,hyper_stack,cline_para,image_para,cline_para.upFactor,flag);
        
        membrane2=coord2image3d(x,y,z,image_para.imsize,CLscale,1);
        mfill=imfill(double(membrane2~=0),'holes');
        CL=[mean(x)' mean(y)' mean(z)'];
        M3=bwdist(~mfill).^2;
        M3=M3/max(M3(:))*CLscale;
        CL = ActiveContourFit(M3, cline_para, CL*CLscale);
        CL=CL/CLscale;
        %%
        if flag.F1==1
            
            xf=x;
            yf=y;
            zf=z;
            
            
            
            F=zeros(row,col,image_para.Fstack_z_size);
            for tt=1:image_para.Fstack_t_size
                % if theres a Fluorescent stack, open it and scale it so that the
                %axis are equal and have its center be the same as the membrane stack
                F1=zeros(row,col,image_para.Fstack_z_size);
                for k = 1:image_para.Fstack_z_size
                    if flag.reverse==0
                        tempImageIdx = image_para.stack_z_size*image_para.stack_t_size+(tt-1)*image_para.Fstack_z_size+k;                       
                    else
                        tempImageIdx = (tt-1)*image_para.Fstack_z_size+k;
                    end
                     tempImage = double(imread([p,'.tif'], 'tif',...
                            tempImageIdx + intensity_correction));
                    % this will be the section that shifts images if there
                    % is metadata
                        if isMetadata
                            filterUsed = totalImageInfo{tempImageIdx}.filterID;
                            % current structure puts the filter number at
                            % the end, such as 641/75 nm (3) is filter 3
                            filterUsed = str2double(filterUsed(end-1)); 
                            F1(:,:,k) = shiftImgFilter(tempImage,filterUsed);
                        else
                            F1(:,:,k) = tempImage;
                        end
                    
%                         % for now, ignore this and continue shifting as
%                         % filter 2
% %                      F1(:,:,k) = shiftImgFilter(tempImage,2);
                end
                
                F1=image_resize(F1,size(F1,1),...
                    size(F1,2),round(image_para.ZF_factor*size(F1,3)));
                Fim_size=size(F1);
                
        pedMaskf=false(size(F1));
        pedMaskf(1:3,:,:)=true;
        pedMaskf(:,1:3,:)=true;
        pedMaskf(end-2:end,:,:)=true;
        pedMaskf(:,end-2:end,:)=true;
        pedMaskf(F1==0)=false;      
        pedistal=median(F1(pedMaskf));
        F1(F1<pedistal)=pedistal;
        F1=normalizeRange(F1);
        if isfield(flag,'F1fix')
        if flag.F1fix
           AAA=convnfft(hyper_stack.*~pedMask,...
           F1(end:-1:1,end:-1:1,end:-1:1).*~pedMaskf);
           [xfix,yfix,zfix]=ind2sub(size(AAA),find(AAA==max(AAA(:))));
           xfix=xfix-ceil(size(AAA,1)/2);
           yfix=yfix-ceil(size(AAA,2)/2);
           zfix=zfix-ceil(size(AAA,3)/2);
           
           
[Xmap,Ymap,Zmap]=ndgrid(1:size(F1,1),1:size(F1,2),1:size(F1,3));
            Xmap=Xmap-xfix;
            Ymap=Ymap-yfix;
            Zmap=Zmap-zfix;
            F1=interp3(F1,Ymap,Xmap,Zmap,'nearest',0);
        end 
        end
        
                
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
                if Fim_size(3)~=image_para.imsize(3)
                    [Ygrid,Xgrid,Zgrid]=meshgrid(1:Fim_size(2),1:Fim_size(1),1:Fim_size(3));
                    Zgrid=Zgrid+(image_para.imsize(3)-Fim_size(3))/2;
                    Ff=interp3(Ygrid,Xgrid,Zgrid,F1,yf,xf,zf,'*linear');
                else
                    Ff=interp3(F1,yf,xf,zf,'*linear');
                end
                [x0,y0,z0,Ff0]=...
                    trimDataEnds(end_slices,xf,yf,zf,Ff);
                I=interp3(hyper_stack,yf, xf,zf,'*linear');
                
                %max radial projection
                
                Ffr=zeros(size(x1));
                r_steps=0.75:.05:1.25;  %radius range
                r_steps=reshape(r_steps,1,1,[]);
                rkx=bsxfun(@minus,xf,CL(:,1)');
                rky=bsxfun(@minus,yf,CL(:,2)');
                rkz=bsxfun(@minus,zf,CL(:,3)');
                rkx=bsxfun(@times,rkx,r_steps);
                rky=bsxfun(@times,rky,r_steps);
                rkz=bsxfun(@times,rkz,r_steps);
                rkx=bsxfun(@plus,rkx,CL(:,1)');
                rky=bsxfun(@plus,rky,CL(:,2)');
                rkz=bsxfun(@plus,rkz,CL(:,3)');
                if Fim_size(3)~=image_para.imsize(3)
                    Ffr=interp3(Ygrid,Xgrid,Zgrid,F1,rky,rkx,rkz);
                    
                else
                    
                    Ffr=interp3(F1,rky,rkx,rkz);
                end
                Ir=interp3(hyper_stack,rky,rkx,rkz);
                Ir=max(Ir,[],3);
                Ffr=max(Ffr,[],3);
                
                end
        end
        
        %% fit polymers
        

if isfield(flag,'polyFit');
if flag.polyFit
%Ff2=interp3(bpass3_jn(F1,.25,[7,7,7]),xf,yf,zf);
polyFit = PolymerFit(F1, xf,yf,zf,Ff,cline_para,image_para );
else
    polyFit=[];
end
else
    polyFit=[];
end     
        %% Calculate curvature and centerline
        %gm and gc are mean and gaussian curvatures, xfyfzf are same as xyz,
        %but closed at the end so that top and bottom values repeat(for
        %visualization of closed cell.
        
        x=x*image_para.nm_per_pixel;
        y=y*image_para.nm_per_pixel;
        z=z*image_para.nm_per_pixel;
        CL=CL*image_para.nm_per_pixel;
        
        
        
        disp(['Finished second round of fits in ',stack_name,'. ',...
            'Elapsed time on this cell: ',num2str(toc)]);
        
        
        start_pt2=start_pt*image_para.nm_per_pixel;
        end_pt2=end_pt*image_para.nm_per_pixel;
        [xf,yf,zf,gc,gm,kplus,kminus]=cellCurvature(x,y,z, Nv,Bv);
        [r,CL,slice_binary,avg_radius]=CLRadiusFix(xf,yf,zf,CL);
        
        [x0,y0,z0,gc0,gm0,kplus0,kminus0,r0]=...
            trimDataEnds(end_slices,xf,yf,zf,gc,gm,kplus,kminus,r);
        
        if flag.F1==1
            if size(Ff)<size(xf,1)
                Ff=cat(1,Ff,Ff(1,:,:));
                Ffr=cat(1,Ffr,Ffr(1,:,:));
            end
            [Ffr0,Ff0,I0]=trimDataEnds(end_slices,Ff,Ffr,I);
        end
        
        %% Calculate surface area, volume and centerline curvature;
        Area=surfArea(xf,yf,zf);
        Volume=surfVolume(xf,yf,zf);
        cellLength=sum(sqrt(sum(diff(CL).^2,2)));
        K=Curvature(CL,5);
        T=Twist(CL,5);
        [pole_left,pole_right] = pole(gm,slice_binary,CL);
        
        
       %% Cluster outputs
       
       cellcoord.x=xf;
       cellcoord.y=yf;
       cellcoord.z=zf;
       cellcoord.CL=CL;
      
       
       cellcurve.kg=gc;
       cellcurve.km=gm;
       cellcurve.k1=kplus;
       cellcurve.k2=kminus;
       cellcurve.twist=T;
       cellcurve.curvature=K;
       
       cellShape.Area=Area;
       cellShape.Volume=Volume;
       cellShape.Length=cellLength;
       cellShape.AvgRadius=avg_radius;
       cellShape.Radius=r;
       cellShape.poleLeft=pole_left;
       cellShape.poleRight=pole_right;
       
       if flag.F1
       cellF.Ff=Ff;
       cellF.Ffr=Ffr;
       cellF.I=I;
       cellF.Ir=Ir;
       else
           cellF=[];
       end
       
       
       
        
        
        
        
        %% Save files

                p2=[p,ptail];
                
                
                save(p2,'p2','cellcoord','cellF','cellcurve','cellShape',...
                    'cline_para','image_para','eg','polyFit')
               
        
        
        %       clear xpe ype zpe xpef ypef zpef xp2ef yp2ef zp2ef xps yps zps xp2sf yp2sf zp2sf xpsf ypsf zpsf
        %       clear xp2 yp2 zp2 xp2e yp2e zp2e xp2s yp2s zp2s Ftest cl centerline
        
        %%
        if Fflag2==1 %%this is for DNA stain or other signal on interior of cell
            error('cellShape:Fflag2','Secondary fluorescent flag set. This setting is currently unavailable');
            
            F=zeros(row,col,image_para.Fstack_z_size);
            
            for k = 1:Fstack_z_size
                F(:,:,k) = double(imread([p,'.tif'], 'tif',image_para.stack_z_size*image_para.stack_t_size*2+k)+intensity_correction);
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
        
        disp(['Finished ',stack_name,'. Elapsed time ',num2str(toc),' seconds.']);
%     catch ME
%         save(fullfile(folderName,[stack_name,'errorCaught.mat']));
%         display(['Error in file ',stack_name,'. Skipping and proceeding with next.']);
%     end
    
end
% end loop over files -----

%% Save the relevent variables


% disp('done');

%% plot
%surf(x0(1:20,90:end),y0(1:20,90:end),z0(1:20,90:end));
%
% surf(xf,yf,zf,gc)
% axis equal
% caxis([min(min(gc))/5,-min(min(gc))/5]);
