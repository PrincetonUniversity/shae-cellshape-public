function [poly,pts] = ActiveContourFitConvMembraneMulti3(hyper_stack, cline_para, cline_initial,xf,yf,zf,Ff,FZ,positions)

% add snake repulsion
% fits snake to image. snake is confined to membrane given by xf yf zf;
%cline initial is given as coords in unwrap pixels
bestcounter=0;
xold=xf;
yold=yf;
zold=zf;
ZX_factor=.6/.8;
upFactor=cline_para.upFactor;
xf=xf*2^(upFactor-1);
yf=yf*2^(upFactor-1);
zf=zf*2^(upFactor-1);
ebest=Inf;
Htemp=normalizeRange(bpass3_jn(hyper_stack,.25,[15,15,15]));
for i=1:max(FZ(:))
Fboundary=bwdist(FZ==i);
[fxb1,fyb1]=gradient(Fboundary);
fxb(:,:,i)=fxb1;
fyb(:,:,i)=fyb1;

end
if upFactor
hyper_stack=interp3(hyper_stack,(upFactor-1));
end
imsize=size(hyper_stack);
options.FourierKernel=true;
options.Power2Flag=false;
%% Load PSF
averageSurfPSF = getappdata(0,'averageSurfPSF');


aPSF=normalizeRange(averageSurfPSF(:,:,:,1));

aPSF=image_resize(aPSF,size(aPSF,1),...
    size(aPSF,2),round(ZX_factor*size(aPSF,3)));
aPSF=normalizeRange(bpass3_jn(aPSF,.25,[15,15,15]));


if upFactor
    aPSF=interpn(aPSF,upFactor-1);
end
%aPSF=padarray(aPSF,[1,1,1]);

% aPSF2=flipdim(aPSF2,1);
% aPSF2=flipdim(aPSF2,2);
%aPSF2=aPSF2/2+flipdim(aPSF2,3)/2;
aPSFfft=fftn(aPSF,size(aPSF)+imsize-[1,1,1]);
% [yPSF,xPSF,zPSF]=gradient(aPSF2,3);
% yPSFfft=fftn(yPSF,size(aPSF2)+imsize-[1,1,1]);
%
% xPSFfft=fftn(xPSF,size(aPSF2)+imsize-[1,1,1]);
% zPSFfft=fftn(zPSF,size(aPSF2)+imsize-[1,1,1]);


% cline_para.stiff=cline_para.stiff*100;
% cline_para.alpha =1*cline_para.stiff;
% cline_para.beta = cline_para.pts*cline_para.stiff;

% initialize snake
rc = cline_initial;
Aout=[];
for iPolymer=1:length(rc);
    %s_interp = 0:0.1:s(end);
    %xyzs = interp1(s, cline_initial, s_interp, 'spline');
    %keep track of where polymers start
    [m n] = size(rc(iPolymer).coord);
    endsIdx(iPolymer)=m;
    % populating the penta diagonal matrix, with the number of polymers
    A = zeros(m,m);
    brow = zeros(1,m);
    brow(1,1:5) = [cline_para.beta, -(cline_para.alpha + 4*cline_para.beta), (2*cline_para.alpha + 6 *cline_para.beta), -(cline_para.alpha + 4*cline_para.beta), cline_para.beta];
    A(1, 1:3) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];
    A(2, 1:4) = [-cline_para.alpha-2*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-4*cline_para.beta, cline_para.beta];
    for i=3:m-2
        A(i,:) = brow;
        brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
    end
    A(m-1, m-3:m) = [cline_para.beta, -cline_para.alpha-4*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-2*cline_para.beta];
    A(m, m-2:m) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];
    Aout=blkdiag(Aout,A);
    
end
endsIdx=cumsum(endsIdx);
startIdx=[1,endsIdx(1:end-1)+1];

[L U] = lu(Aout + cline_para.gamma .* eye(length(Aout)));
%Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization

% find the center line of the cell

% two pass algorithm

% vidObj = VideoWriter('tmp.avi', 'Motion Jpeg AVI');
% vidObj.FrameRate = 15;
% vidObj.Quality = 95;
% open(vidObj);
%k: stack_z_size index
%obtain stack_z_size(t)

%% filter image to "normalize" intensity
im=hyper_stack;
%gradient



%% prepare for cline_para.iterations
rc_avg = zeros(size(rc));
counter = 0;
Kscale=1;
Kr=1;



%calculate u and v tangent vectors
[mf,nf]=size(xf);
l=3;
x=[xf;xf];
y=[yf;yf];
z=[zf;zf];
[Xu(:,:,1),Xv(:,:,1)]=gradient(x,l);
[Xu(:,:,2),Xv(:,:,2)]=gradient(y,l);
[Xu(:,:,3),Xv(:,:,3)]=gradient(z,l);

Xu(1:2*l,:,:)=Xu(mf+1:mf+2*l,:,:);
Xv(1:2*l,:,:)=Xv(mf+1:mf+2*l,:,:);
Xu=Xu(1:mf,:,:);
Xv=Xv(1:mf,:,:);




%% moving the snake in each iteration
for i=1:cline_para.iterations;
    %calculate image gradient force on the nodes
    
    
    if 1                % actual interpolation, slow!!!
        wMask=ones(size(hyper_stack));
        CLim=zeros(size(hyper_stack));
        xyzs=[];
        rcall=[];
        uvec=[];
        vvec=[];
        totLength=0;
        Npolymer=length(rc);
        
        for iPolymer=1:Npolymer
            Ccoord=rc(iPolymer).coord(:,1);
            Rcoord=rc(iPolymer).coord(:,2);
            polymer3(iPolymer).xs=interp2(xf,Ccoord,Rcoord,'*spline');
            polymer3(iPolymer).ys=interp2(yf,Ccoord,Rcoord,'*spline');
            polymer3(iPolymer).zs=interp2(zf,Ccoord,Rcoord,'*spline');
            polymer3(iPolymer).length=calculate_length(polymer3(iPolymer).xs,...
                polymer3(iPolymer).ys,polymer3(iPolymer).zs);
            
            
            
            uv=[interp2(Xu(:,:,1),Ccoord,Rcoord,'*spline'),...
                interp2(Xu(:,:,2),Ccoord,Rcoord,'*spline'),...
                interp2(Xu(:,:,3),Ccoord,Rcoord,'*spline')];
            vv=[interp2(Xv(:,:,1),Ccoord,Rcoord,'*spline'),...
                interp2(Xv(:,:,2),Ccoord,Rcoord,'*spline'),...
                interp2(Xv(:,:,3),Ccoord,Rcoord,'*spline')];
            polymer3(iPolymer).uvec=normr(uv);
            polymer3(iPolymer).vvec=normr(vv);
            
            C1=coord2image3d(polymer3(iPolymer).xs,polymer3(iPolymer).ys,...
                polymer3(iPolymer).zs,size(hyper_stack),1,1);
            
            convFit1=normalizeRange(convnfft(C1, aPSFfft,'same',1:3,options));
            
   %         wMask(convFit1>.5 & hyper_stack>median(Ff(FZ==iPolymer)))=0;
            convFit1=convFit1*quantile(Ff(FZ==iPolymer),.95);
            
            CLim=CLim+convFit1;
            xyzs=cat(1,xyzs,[polymer3(iPolymer).xs,polymer3(iPolymer).ys,polymer3(iPolymer).zs]);
            rcall=cat(1,rcall,rc(iPolymer).coord);
            uvec=cat(1,uvec,polymer3(iPolymer).uvec);
            vvec=cat(1,vvec,polymer3(iPolymer).vvec);
            
            
            polymer3(iPolymer).length=sum(sqrt(sum(diff([polymer3(iPolymer).xs,polymer3(iPolymer).ys,polymer3(iPolymer).zs],[],1).^2,2)));
            totLength=totLength+polymer3(iPolymer).length;
        end
     convFit=CLim;
        if i<300
            deltaIm=convFit-Htemp;
        else
            deltaIm=convFit-hyper_stack;
            
        end
        deltaIm=deltaIm.*wMask;
        deltaIm(hyper_stack==0)=0;
        
        e=sum(deltaIm(:).^2);
        if e<ebest && i>15;
            ebest=e;
            rc_best=rc;
            bestcounter=0;            
            if isfield(cline_para,'thermal')
            %    cline_para.kappa=cline_para.kappa*1.02;
                
                cline_para.thermal=cline_para.thermal*.999;
            end
            
        end
        bestcounter=bestcounter+1;
        if bestcounter>50
            display(['Limit reached after ' num2str(i) ' iterations']);
            break
        end
        
        
        if i>15
          deltaImConv=convnfft(deltaIm,aPSFfft,'same',1:3,options)/20;
            
            [fy,fx,fz]=gradient(deltaImConv*1);
            cline_para.stretch_ends_flag=1;
        else
            [fy,fx,fz]=gradient(deltaIm*1);
            cline_para.stretch_ends_flag=1;
            deltaImConv=deltaIm;
            
        end
        fs = [interp3(fx,xyzs(:,2),xyzs(:,1), xyzs(:,3),'*spline'),...
            interp3(fy,xyzs(:,2),xyzs(:,1), xyzs(:,3),'*spline'), ...
            interp3(fz,xyzs(:,2),xyzs(:,1), xyzs(:,3),'*spline')];
        
        fs=fs*cline_para.kappa;
        fs(isnan(fs)) = 0;
    else
        s_ind = round(xyzs);
        fs = zeros(m, 3);
        in_image_flag = all((s_ind>0)&(bsxfun(@le, s_ind, [col, row, stack_z_size])), 2);
        gradient_ind = sub2ind([row, col, stack_z_size], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1), s_ind(in_image_flag, 3));
        fs(in_image_flag, :) = -[fx(gradient_ind), fy(gradient_ind), fz(gradient_ind)];
    end
    
    %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
    % viscos force between frames on all nodes:
    
    % stretch or shrink the snake at the ends
    if cline_para.stretch_ends_flag
        for iPolymer=1:length(rc);
            
            s_head = xyzs(startIdx(iPolymer),:) - xyzs(startIdx(iPolymer)+4, :);
            s_head = s_head/norm(s_head);
            Ftipstart=interp3(deltaImConv,xyzs(startIdx(iPolymer),2),...
                xyzs(startIdx(iPolymer),1),xyzs(startIdx(iPolymer),3),'*linear');
            
            s_tail = xyzs(endsIdx(iPolymer), :) - xyzs(endsIdx(iPolymer) - 3, :);
            s_tail = s_tail/norm(s_tail);
            
            Ftipend=interp3(deltaImConv,xyzs(endsIdx(iPolymer),2),...
                xyzs(endsIdx(iPolymer),1),xyzs(endsIdx(iPolymer),3),'*linear');
            
            %only if very straigth
            Fends(1) = -Ftipstart*cline_para.stretching_force_factor;
            Fends(2) =Ftipend*cline_para.stretching_force_factor;
            if polymer3(iPolymer).length<.6
                Fends(1)=Fends(1)-4;
                Fends(2)=Fends(2)+4;
            end
            
            fstretch=interp1(0:1,Fends,linspace(0,1,endsIdx(iPolymer)-startIdx(iPolymer)+1));
            [~,tVector]=gradient(xyzs(startIdx(iPolymer):endsIdx(iPolymer),:),3);
            tVector=normr(tVector);
            fstretch=bsxfun(@times,tVector,fstretch');
            fs(startIdx(iPolymer):endsIdx(iPolymer),:)=fs(startIdx(iPolymer):endsIdx(iPolymer),:)+fstretch;
        end
    else
        Ftipstart=0;
        Ftipend=0;
        % xyzstart=xyzs(1,:);
        % xyzend=xyzs(end,:);
    end
    
    if i>Inf
        Frepulsion=zeros(size(xyzs));
        for ipoly=2:Npolymer
            for jpoly=1:ipoly-1
                xmat=bsxfun(@minus,xyzs(startIdx(ipoly):endsIdx(ipoly),1)',xyzs(startIdx(jpoly):endsIdx(jpoly),1));
                ymat=bsxfun(@minus,xyzs(startIdx(ipoly):endsIdx(ipoly),2)',xyzs(startIdx(jpoly):endsIdx(jpoly),2));
                zmat=bsxfun(@minus,xyzs(startIdx(ipoly):endsIdx(ipoly),2)',xyzs(startIdx(jpoly):endsIdx(jpoly),2));
                
                dmat=sqrt(xmat.^2+ymat.^2+zmat.^2);
                xmat=xmat./dmat;
                ymat=ymat./dmat;
                zmat=zmat./dmat;
                
                Fmat=(1-dmat/cline_para.objsize).^2;
                Fmat(dmat>cline_para.objsize)=0;
                
                Frepulsionx=Fmat.*xmat;
                Frepulsiony=Fmat.*ymat;
                Frepulsionz=Fmat.*zmat;
                Frepulsion(startIdx(ipoly):endsIdx(ipoly),:)=...
                    Frepulsion(startIdx(ipoly):endsIdx(ipoly),:)-...
                    [sum(Frepulsionx,1)',sum(Frepulsiony,1)',sum(Frepulsionz,1)'];
                
                Frepulsion(startIdx(jpoly):endsIdx(jpoly),:)=...
                    Frepulsion(startIdx(jpoly):endsIdx(jpoly),:)+...
                    [sum(Frepulsionx,2),sum(Frepulsiony,2),sum(Frepulsionz,2)];
                
                
                
            end
        end
        fs=fs+Frepulsion*20;
    end
    
    fu=dot(fs,uvec,2);
    fv=dot(fs,vvec,2);
    fs=[fu,fv];
    
    
    if isfield(cline_para,'thermal')
        
        for iPolymer=1:length(rc)
            thermal=max(max(abs(fs(startIdx(iPolymer):endsIdx(iPolymer),:))))*cline_para.thermal;
            
%                         Ftrans=thermal*(rand(1,2)-[.5,.5])/5;
%                         fs(startIdx(iPolymer):endsIdx(iPolymer),:)=...
%                             bsxfun(@plus,fs(startIdx(iPolymer):endsIdx(iPolymer),:),Ftrans);
%             
            Frot=normalizeRange((startIdx(iPolymer):endsIdx(iPolymer)));
            Frot=Frot-mean(Frot);
            v_perp=rc(iPolymer).coord(1,:)-rc(iPolymer).coord(end,:);
            v_perp=normr(v_perp);
            v_perp=[-v_perp(2),v_perp(1)];
            Frot=Frot'*v_perp*thermal*(rand(1)-.5)*5;
            fs(startIdx(iPolymer):endsIdx(iPolymer),:)=...
                fs(startIdx(iPolymer):endsIdx(iPolymer),:)+Frot;
            
        end
        
    end
    %calculate worm confinement force
    
    if any(~FZ(:))
        fsb=zeros(size(fs));
        for iPolymer=1:max(FZ(:))            
        fsb(startIdx(iPolymer):endsIdx(iPolymer),:) = [interp2(fxb(:,:,iPolymer)...
            ,rcall(startIdx(iPolymer):endsIdx(iPolymer),1),...
            rcall(startIdx(iPolymer):endsIdx(iPolymer),2),'*spline'),...
            interp2(fyb(:,:,iPolymer),rcall(startIdx(iPolymer):endsIdx(iPolymer),1),...
            rcall(startIdx(iPolymer):endsIdx(iPolymer),2),'*spline')];
       
        end
         fs=fs+fsb*50;
        
    end
    
    
    
    
    
    %calculate the new position of snake
    new_rc =U\(L\(cline_para.gamma*rcall- fs));
    delta_rc=new_rc-rcall;
    rcall=rcall+delta_rc;
    for iPolymer=1:length(rc)
        rc(iPolymer).coord=respace(rcall(startIdx(iPolymer):endsIdx(iPolymer),:));
    end
    
    
    if cline_para.show_flag==1
            [i,e, totLength*80/5,ebest]

  %      memIm=double(interp3(Htemp,yf,xf,zf));
        % if nargin>7
        %         cell=surf(xf,yf,zf,Ff,'edgecolor','none','facealpha',.2);
        % else
        %     Ff=interp3(hyper_stack,yf,xf,zf);
        %             cell=surf(xf,yf,zf,Ff,'edgecolor','none','facealpha',.2);
        %
        % end
        %
        subplot(1,2,1);
        cell=surf(xf,yf,zf,Ff,'facelighting','phong','facecolor','inter','edgecolor','none','facealpha',1);

        %         %      surf(zf,xf,yf,memIm,'edgecolor','none');
        
        hold on
        if nargin>8
            plot3(positions(:,1),positions(:,2),positions(:,3),'xg','linewidth',3)
        end
        % colorbar
        for iPolymer=1:max(FZ(:))            
               pts=plot3(xyzs(startIdx(iPolymer):endsIdx(iPolymer),1)...
                   ,xyzs(startIdx(iPolymer):endsIdx(iPolymer),2)...
                   ,xyzs(startIdx(iPolymer):endsIdx(iPolymer),3),'black','linewidth',3);

        end
        
        
        %    scatter3(xyzs(:,3),xyzs(:,1),xyzs(:,2));
        
        xlim([0,50])
        %       plot(zf(:,20),yf(:,20));
        %
        %               scatter(xyzs(:,3),xyzs(:,2));
        %       hold on
    %    fs2=bsxfun(@times,fs(:,1),uvec)+bsxfun(@times,fs(:,2),vvec);
       % quiver3(xyzs(:,1),xyzs(:,2),xyzs(:,3),fs2(:,1),fs2(:,2),fs2(:,3))
        %  quiver3(xyzs(:,1),xyzs(:,2),xyzs(:,3),delta_xyzs(:,1),delta_xyzs(:,2),delta_xyzs(:,3))
        hold off
        axis equal
        axis off
                subplot(1,2,2);
                set(gcf,'color','w');

        imagesc(Ff);
        hold on
        for iPolymer=1:max(FZ(:))            

        plot(rcall(startIdx(iPolymer):endsIdx(iPolymer),1),...
            rcall(startIdx(iPolymer):endsIdx(iPolymer),2),'black','linewidth',3);

         plot(rcall(startIdx(iPolymer):endsIdx(iPolymer),1),...
            rcall(startIdx(iPolymer):endsIdx(iPolymer),2)-size(Ff,1)/3,'black','linewidth',3);
       
                 plot(rcall(startIdx(iPolymer):endsIdx(iPolymer),1),...
            rcall(startIdx(iPolymer):endsIdx(iPolymer),2)-2*size(Ff,1)/3,'black','linewidth',3);
        end
        
        hold off
        set(gca,'ydir','normal')
        set(gca,'ytick',linspace(1,size(Ff,1)/3,3));
        set(gca,'yticklabel',({'0', 'pi', '2pi'}));
        set(gca,'xticklabel',[]);
        ylim([1,size(Ff,1)/3]);
        %imagesc(max(deltaIm,[],3));
        pause(.1)
    end
end

rc=(rc_best);
xyzs=[];
rcall=[];
for iPolymer=1:length(rc)
    Ccoord=rc(iPolymer).coord(:,1);
    Rcoord=rc(iPolymer).coord(:,2);
    xs=interp2(xf,Ccoord,Rcoord,'*spline');
    ys=interp2(yf,Ccoord,Rcoord,'*spline');
    zs=interp2(zf,Ccoord,Rcoord,'*spline');
    poly(iPolymer).xs=xs;
    poly(iPolymer).ys=ys;
    poly(iPolymer).zs=zs;
    poly(iPolymer).rc=rc(iPolymer).coord;
    poly(iPolymer).length=calculate_length(xs,ys,zs);
    poly(iPolymer).Curvature=Curvature([xs,ys,zs],5);
    poly(iPolymer).Twist=Twist([xs,ys,zs],5)';
    poly(iPolymer).AngleAvg=mean(atan2(poly(iPolymer).Twist,poly(iPolymer).Curvature));
    
    xyzs=cat(1,xyzs,[xs,ys,zs]);
    
    
end
pts=xyzs;


end


function rc=respace(rc)
L=length(rc)-1;
s=[0,cumsum(sqrt(diff(rc(:,1)).^2+diff(rc(:,2)).^2))']';

rc(:,1)=interp1(s,rc(:,1),0:max(s)/L:max(s),'spline');
rc(:,2)=interp1(s,rc(:,2),0:max(s)/L:max(s),'spline');
end

function L=calculate_length(x,y,z)
dx=diff(x);dy=diff(y);dz=diff(z);
dl=sqrt(dx.^2+dy.^2+dz.^2);
L=sum(dl);
end


% close(vidObj);

% save data:
% save('xyzs_all.mat', 'xyzs_all');
% clearvars -except xyzs_all hyper_stack stack_z_size_max stack_t_size;