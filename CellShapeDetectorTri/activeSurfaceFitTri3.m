
function [F,eg]=activeSurfaceFitTri3(F,hyper_stack,cline_para,image_para,flag)
VOXELS_COUNT_A = 15;




% writerObj = VideoWriter('cell4HardtofitTri');
% open(writerObj);

Wimage=1;%1./sqrt(hyper_stack);
Wimage(Wimage>50)=50;
H3=hyper_stack>.1;
vertices_avg=0;
row=image_para.imsize(1);
col=image_para.imsize(2);
stack_z_size=image_para.imsize(3);
imsize=image_para.imsize;
ZX_factor=image_para.ZX_factor;
options.FourierKernel=true;
options.Power2Flag=false;
zeroMask=(hyper_stack>=0.15&hyper_stack<=.5);


%% Load PSF
averageSurfPSF = getappdata(0,'averageSurfPSF');
if ~isfield(cline_para,'upFactor')
    upFactor=cline_para.upFactor;
else
    upFactor=2;
end

if flag.psf4d==1
    aPSF2=normalizeRange(averageSurfPSF);
    aPSF=aPSF2.^(1/upFactor);
    for i=1:size(aPSF,4)
        aPSF(:,:,:,i)=image_resize(aPSF2(:,:,:,i),size(aPSF2,1),...
            size(aPSF2,2),round(ZX_factor*size(aPSF2,3)));
    end
else
%    %%%%%%%%%
%     aPSF=normalizeRange(averageSurfPSF(:,:,:,1)-.01);
%     %%%%%%%%%%
aPSF=averageSurfPSF(:,:,:,1);
    aPSF=image_resize(aPSF,size(aPSF,1),...
        size(aPSF,2),round(ZX_factor*size(aPSF,3)));
end

aPSF2=aPSF;
aPSF=image_resize(aPSF,size(aPSF,1)*upFactor,size(aPSF,2)*upFactor,...
    size(aPSF,3)*upFactor);
%aPSF=padarray(aPSF,[1,1,1]);

% aPSF2=flipdim(aPSF2,1);
% aPSF2=flipdim(aPSF2,2);
% aPSF2=aPSF2/2+flipdim(aPSF2,3)/2;
aPSFfft=fftn(aPSF,size(aPSF)+imsize*upFactor-[1,1,1]);

[yPSF,xPSF,zPSF]=gradient(aPSF2,1);
aPSFfft2=fftn(aPSF2,size(aPSF2)+imsize-[1,1,1]);
yPSFfft=fftn(yPSF,size(aPSF2)+imsize-[1,1,1]);

xPSFfft=fftn(xPSF,size(aPSF2)+imsize-[1,1,1]);
zPSFfft=fftn(zPSF,size(aPSF2)+imsize-[1,1,1]);




%% Create initial membrane outline and simulated image


% end_slices=max(end_slices,start_slices);
% r_vector=cellShellVectors(x1,y1,z1,end_slices);
Nvec=vertnorm_dw(F.faces,F.vertices);

[~,vertAreaI,~]=triArea(F);

[F2,vertAreaI2]=triInterp2(F,vertAreaI,2);

membrane=coord2image3d(F2...
    ,imsize,upFactor,1,vertAreaI2);
vertArea=vertAreaI2;
if cline_para.gradflag
    membrane=fillMembrane(membrane,F2);
end


% vertices=[reshape(x1,numel(x1),1),reshape(y1,numel(y1),1),...
%     reshape(z1,numel(z1),1)];
%vertices=vertices;

if flag.fftconv==1
    convFit=convnfft(membrane, aPSFfft,'same',1:3,options);
    
    
else
    convFit=conv3sparse(membrane, aPSF,'same');
end

%% 3D Forces
%Calculatetr

%FL and FU are te upper and lower triangular factorizations that can
%be inverted
[FL,FU]=Finternal3Dmat(F,cline_para);
[FLz,FUz]=Finternal3Dzmat(F,cline_para);



%% Cut off points that lie outside the image, scale simulated image
% to be similar to actual image
Xpts=(1:upFactor:upFactor*row)+1;
Ypts=(1:upFactor:upFactor*col)+1;
Zpts=(1:upFactor:upFactor*stack_z_size)+1;
[Ypts,Xpts,Zpts]=meshgrid(Ypts,Xpts,Zpts);

Xpts2=(0:1/upFactor:row)+1/upFactor/2;
Ypts2=(0:1/upFactor:col)+1/upFactor/2;
Zpts2=(0:1/upFactor:stack_z_size)+1/upFactor/2;
[Ypts2,Xpts2,Zpts2]=meshgrid(Ypts2,Xpts2,Zpts2);

%%

convFitSub=interp3(convFit,Ypts,Xpts,Zpts);
vertices=F.vertices;
%vertices(all(isnan(vertices),2),:)=[];
vertices(vertices(:,1)>row,1)=row;
vertices(vertices(:,2)>col,2)=col;
vertices(vertices(:,3)>stack_z_size,3)=stack_z_size;
vertices(vertices(:,1)<1,1)=1;
vertices(vertices(:,2)<1,2)=1;
vertices(vertices(:,3)<1,3)=1;


vtemp=vertices;
vtemp(all(isnan(vtemp),2),:)=[];
imAvg=mean(hyper_stack(sub2ind(size(hyper_stack),round(vtemp(:,1)),...
    round(vtemp(:,2)),round(vtemp(:,3)))));
fitAvg=mean(convFitSub(sub2ind(size(hyper_stack),round(vtemp(:,1)),...
    round(vtemp(:,2)),round(vtemp(:,3)))));
K=mean(convFitSub(H3).*hyper_stack(H3))/mean(convFitSub(H3).^2);
if isnan(K)
    K=1;
end

Kr=1;


deltaMembrane=zeros(size(membrane));
membrane=membrane*K;
convFitSub=convFitSub*K;
convFit=convFit*K;
%% Update fits
best=Inf;
counter=0;
for i=1:cline_para.iterations3d
    
    membrane=membrane+deltaMembrane;
    if flag.fftconv==1
        
        convFit=convnfft(membrane, aPSFfft,'same',1:3,options);
        
    else
        deltaConv=conv3sparse(deltaMembrane,aPSF,'same');
        convFit=convFit+deltaConv;
    end
    convFitSub=interp3(convFit,Ypts,Xpts,Zpts);
    
        K=mean(convFitSub(H3).*hyper_stack(H3))/mean(convFitSub(H3).^2)*.8;

    convFitSub=convFitSub*K;
    if cline_para.gradflag
      %  scale bright portions of both images to be more similar
      maxMask=hyper_stack>.8 & convFitSub>.5;
      approxObjSize=max(max(max(bwdist(~membrane))));
      if nnz(maxMask)<VOXELS_COUNT_A
         warning('activeSurfaceFit:voxelCountTooLow','Too few voxels above threshold for normalization');
         maxMask = hyper_stack>0.6 & convFitSub>0.4;
         if nnz(maxMask)<VOXELS_COUNT_A
         warning('activeSurfaceFit:voxelCountTooLow2','Too few voxels above secondary threshold for normalization');
         maxMask = hyper_stack>0.35 & convFitSub>0.2;
         end
      end
        maxMask=imdilate(maxMask,true(3,3,3));
        maxMean=mean(convFitSub(maxMask));
        if approxObjSize>8
        convFitSub=(convFitSub.*(mean(hyper_stack(maxMask))))/maxMean*1;
        else
        convFitSub=(convFitSub.*(mean(hyper_stack(maxMask))))/maxMean*.8;
        end
    end
    convFitSub=(convFitSub);
    deltaIm=(convFitSub-hyper_stack);
   % deltaIm(abs(deltaIm)<.01)=0;
    if flag.fftconv==0
        if i==1
            deltaImE=conv3sparse(deltaIm.*(convFitSub~=0),aPSF2,'same');
            deltaImOld=zeros(size(deltaImE));
        end
        deltaImD=deltaIm-deltaImOld;
        deltaImE=deltaImE+ conv3sparse(deltaImD,aPSF2,'same');%./areaCorrection;
        deltaImOld=deltaIm;
        
    end
    
    
    if i<15 && ~cline_para.gradflag
        deltaImE=deltaIm;
        %   [Fy,Fx,Fz]=gradient(convnfft(deltaIm/10,aPSFfft,'same',1:3,options),3);
        [Fy,Fx,Fz]=gradient(-hyper_stack);
        
        fs=[interp3(Fx,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
            interp3(Fy,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
            interp3(Fz,vertices(:,2),vertices(:,1),vertices(:,3),'*linear')];
    else
        
        %         Fz=convnfft(deltaIm,zPSFfft,'same',1:3,options)/20;
        %         Fx=convnfft(deltaIm,xPSFfft,'same',1:3,options)/20;
        %         Fy=convnfft(deltaIm,yPSFfft,'same',1:3,options)/20;
        if cline_para.gradflag==0
            deltaImConv=convnfft(Wimage.*deltaIm,aPSFfft2,'same',1:3,options);
            
            [Fy,Fx,Fz]=(gradient(deltaImConv));
            fs=[interp3(Fx,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
                interp3(Fy,vertices(:,2),vertices(:,1),vertices(:,3),'*linear'),...
                interp3(Fz,vertices(:,2),vertices(:,1),vertices(:,3),'*linear')];
            if mod(i,6)-1==0 || ~exist('Nvec','var')
                Nvec=vertnorm_dw(F.faces,F.vertices);
            end
            
%             if i==60
%                 cline_para.kappa3d=1.5*cline_para.kappa3d;
%             end
        else
        %    deltaIm=deltaIm.*zeroMask;
            
            deltaImConv=convnfft(Wimage.*deltaIm,aPSFfft2,'same',1:3,options);
            %deltaImConv=deltaIm*100;
            fs=interp3(deltaImConv,vertices(:,2),vertices(:,1),vertices(:,3),'*linear');
            if mod(i,6)-1==0 || ~exist('Nvec','var')
                Nvec=vertnorm_dw(F.faces,F.vertices);
            end
            fs=bsxfun(@times,Nvec,fs);
            
            
        end
    end
    
    if i==11
        eg=convFitSub;
        Fbest=F;
        %         bestx=x;
        %         besty=y;
        %         bestz=z;
        counter=0;
        % reset the energy to infinity after initial relaxation
        best = inf;
    end
    %     dK=interp3(deltaIm,Ypts2,Xpts2,Zpts2);
    %     dK=dK(boolean(imdilate(membrane>0,true(3,3,3))));
    dK=deltaIm;
    %(nansum(dK(:).^2))
%     % the best energy must a strictly positive number. Correct this if it
%     % is ever not true
%     if isnan(best) || best<=0
%         best = inf;
%     end
    if (nansum(dK(:).^2)<best&& i>11) % || i==11
        eg=convFitSub;
        Fbest=F;
        %         bestx=x;
        %         besty=y;
        %         bestz=z;
        counter=0;
        % calculate the image energy as the nansum or sum?
        % BPB changes to nansum 20140506
        best=nansum(dK(:).^2);
        
    elseif i>11
        counter=counter+1;
    end
    if counter==25;
        break
    end
    
    
    %calculate forces at each point by interpolating
    
    vertices(isnan(vertices))=0;
    fs(isnan(fs))=0;
    
    %Find distance to move by matrix division
    
    
    if i<10 && cline_para.gradflag
        new_vertices(:,1:2) = FU\(FL \...
            (cline_para.gamma*vertices(:,1:2)));
        new_vertices(:,3) = FUz\(FLz \ ...
            (cline_para.gamma*vertices(:,3)));
        delta_vert=new_vertices-vertices;
    else
        new_vertices(:,1:2) = FU\(FL \...
            (cline_para.gamma*vertices(:,1:2)-cline_para.kappa3d*fs(:,1:2)));
        new_vertices(:,3) = FUz\(FLz \ ...
            (cline_para.gamma*vertices(:,3)-cline_para.kappa3d*fs(:,3)));
        delta_vert=new_vertices-vertices;
        
    end
    
    
    %  delta_vert(abs(delta_vert)>1)=0;
    %     delta_vert(1:200,:)=0;
    %         delta_vert(end-200:end,:)=0;
    delta_vert(vertices==0)=0;
    if ~cline_para.gradflag% && i<10
        if i>0
            
            %find change in areas
            F2=F;
            F2.vertices=new_vertices;

            [~,vertArea2]=triInterp2(F2,[],2);
            [~,vertAreaI2]=triInterp2(F,[],2);
            
            deltaA=vertArea2-vertAreaI2;
            
            if i<10
                delta_vert=zeros(size(F.vertices));%bsxfun(@times,dot(Nvec,delta_vert,2),Nvec);
            elseif mod(i,10 )~=0;
                delta_vert=bsxfun(@times,dot(Nvec,delta_vert,2),Nvec);
            end
        end
    else

    end
    new_vertices=vertices+delta_vert;

    vertices=new_vertices;
  %  Average over the last 20 iterations

    vtemp=vertices;
    vtemp(all(isnan(vtemp),2),:)=[];
    vtemp(all(vtemp==0,2),:)=[];
    F.vertices=vertices;
    F.vertices(F.vertices==0)=NaN;
    F2=triInterp2(F,[],2);


    
    if i>10 && ~cline_para.gradflag% && ~mod(i,10)
        deltaA(deltaA>0)=0;
       vertAreaI2=vertAreaI2-deltaA;
%        [~,vertAreaI2,~]=triArea(F2);
        membrane2=coord2image3d(F2...
            ,imsize,upFactor,1,vertAreaI2);
    elseif ~cline_para.gradflag
        membrane2=coord2image3d(F2...
            ,imsize,upFactor,1,vertAreaI2);
    end
    
    if cline_para.gradflag
        
        F3=triInterp2(F2);
        membrane2=coord2image3d(F3...
            ,imsize,upFactor,1);
        membrane2=fillMembrane(membrane2,F2);
        
    end
    
    
    
    if cline_para.show_flag>0
        % display some interesting numbers
        [(nansum(dK(:).^2)),i,best]
        % seems to have problems making a figure if one doesn't already
        % exist. 
        % figure(gcf);
        if exist('h','var') && all(ishandle(h)), delete(h); end
        quiver3(F.vertices(:,1),F.vertices(:,2),F.vertices(:,3)...
            ,delta_vert(:,1),delta_vert(:,2),delta_vert(:,3));
        
        hold on
        quiver3(F.vertices(:,1),F.vertices(:,2),F.vertices(:,3)...
            ,fs(:,1),fs(:,2),fs(:,3),'r');
        
        
        h = patch(F,'facecolor','w','facealpha',.8);  axis equal;  view(3);
        %   h=patch(F2,'facevertexCdata',vertAreaI,'facecolor','interp','facealpha',.5,'edgecolor','none');  axis equal;  view(3);
        % xlim([-7,50]);ylim([15,60]);zlim([5,35]);
        drawnow;
        hold off
        % frame = getframe;
        % writeVideo(writerObj,frame);
        
        if cline_para.show_flag == 3
            % this will save the current surface to the base workspace for
            % later usage
            if i == 1
                % clear and save the first set
                clear savedSurface;
                savedSurface.fs = fs;
                savedSurface.F = F;
                savedSurface.delta_vert = delta_vert;
                savedSurface.imageEnergy = nansum(dK(:).^2);
                % save it in the base workspace
                assignin('base','savedSurface',savedSurface);
                
            else
                % retrieve original save
                savedSurface = evalin('base','savedSurface');
                % append current surface
                savedSurface(i).fs = fs;
                savedSurface(i).F = F;
                savedSurface(i).delta_vert = delta_vert;
                savedSurface(i).imageEnergy = nansum(dK(:).^2);
                % save
                assignin('base','savedSurface',savedSurface);
                
            end
        end
        
    end
    
    membrane2=membrane2;%.*K;
    % rarely there is an issue with a sizing mismatch, postpend a bunch of
    % zeros to make there be a match
    if not(isequal(size(membrane2),size(membrane)))
       if numel(membrane2)>numel(membrane)
            membrane = padarray(membrane,[10,10,10],0,'post');
            
            membrane = membrane(1:size(membrane2,1),:,:);
            membrane = membrane(:,1:size(membrane2,2),:);
            membrane = membrane(:,:,1:size(membrane2,3));
       else
           membrane2 = padarray(membrane2,[10,10,10],0,'post');
           
            membrane2 = membrane2(1:size(membrane,1),:,:);
            membrane2 = membrane2(:,1:size(membrane,2),:);
            membrane2 = membrane2(:,:,1:size(membrane,3));
       end
    end
       % check to see if there is a size mismatch
    assert(numel(membrane2)==numel(membrane));
    deltaMembrane=membrane2-membrane;
    [~,Elength]=EdgeLength(F);
    % resurface sample every 10 iterations
    if  (mod(i,10)==0 && i<cline_para.iterations3d-21) ||...
            max(Elength)>cline_para.upperLimit+1
        [F,vertAreaI]=surfaceResample4(F,cline_para.upperLimit,cline_para.lowerLimit,vertAreaI);
        
        [FL,FU]=Finternal3Dmat(F,cline_para);
        [FLz,FUz]=Finternal3Dzmat(F,cline_para);
        Nvec=vertnorm_dw(F.faces,F.vertices);
        vertices=F.vertices;
        new_vertices=vertices;
    end
    %     imagesc(squeeze(membrane(:,42,:))');axis equal;
    %     pause(.1)
end
F=Fbest;% x=bestx;
% y=besty;
% z=bestz;



% [kx,ky,kz]=meshgrid(1:201,1:201,1:201);
% kx=kx-101;
% ky=ky-101;
% kz=kz-101;
%K=(kx.^2+ky.^2+kz.^2);
%E=convnfft(membrane2,K,'same');
%E=normalizeRange(E);
%E=E.*mfill;
%E=E-min(E(E~=0));
%E=E+~mfill*1;
%E=E*100;
F=faceNormals(F);




