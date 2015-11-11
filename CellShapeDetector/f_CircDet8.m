function [xys_all,c1,c2]= f_CircDet8(stack, cline_para)
%
% uses active contour to find circular objects in a series of images
%representing a slices of a cylindrical object.
%stack = V;
%%
minR=5;
%stack=permute(stack,[3,1,2]);
if cline_para.movie_flag==1
    writerObj = VideoWriter(cline_para.movie_title);
    open(writerObj);
end
deltacount=[0,0];

if size(size(stack))==3
    stack1=zeros(size(stack),1);
    stack1(:,:,1)=stack;
    stack=stack1;
end


[row, col, stack_z_size] = size(stack);
kernal=fspecial('sobel');
time_idx = 1:stack_z_size;
lowResidue=Inf;
Rcounter=0;
% initialize snake
if ~isfield(cline_para,'initial')
theta_interp = -pi+2*pi/cline_para.pts:2*pi/cline_para.pts:pi;
n_vector=zeros(cline_para.pts,2);
n_vector(:,1)=cos(theta_interp);
n_vector(:,2)=sin(theta_interp);

xys(:,1) = cline_para.center(1)+cline_para.radius*cos(theta_interp);
xys(:,2)= cline_para.center(2)+cline_para.radius*sin(theta_interp);
end

[m n] = size(xys);

% populating the penta diagonal matrix
A = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.beta, -(cline_para.alpha + 4*cline_para.beta), (2*cline_para.alpha + 6 *cline_para.beta), -(cline_para.alpha + 4*cline_para.beta), cline_para.beta];
brow=circshift(brow',-2)';
for i=1:m
    A(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end


[L U] = lu(A + cline_para.gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization

% find the center line of the cell
xys_all = zeros([size(xys), length(time_idx)]);
xys_predicted = xys_all;
im = zeros(row, col);
% two pass algorithm

% vidObj = VideoWriter('tmp.avi', 'Motion Jpeg AVI');
% vidObj.FrameRate = 15;
% vidObj.Quality = 95;
% open(vidObj);

%%

for pass = 1
    %k: stack_z_size index
    for k = 1:length(time_idx)
        %%
        %obtain stack_z_size(t)
        
        if pass==2
            k=length(time_idx)-k+1;
            cline_para.group_iterations=1;
            xys=xys_predicted(:,:,k);
        end
        
        im = stack(:,:,k);
        %% filter image to "normalize" intensity
        
        
        %         if cline_para.gradflag==1;
        %             imMask=(im~=0);
        %             imMask=smooth2a(imMask,2,2);
        %             imMask=(imMask==1);
        %             hoIm=conv2(im,kernal,'same');
        %             vertIm=conv2(im,kernal','same');
        %             im=sqrt(hoIm.^2+vertIm.^2).*(im-min(min(im)))./(max(max(im))-min(min(im)));
        %             im=imMask.*im;
        %
        %
        %         im(1:2,:)=0;
        %         im(:,1:2)=0;
        %         im(row-1:row,:)=0;
        %         im(:,col-1:col)=0;
        %         end
        
        %             % light blur
        x = -12:12;
        h = exp(-x.^2/2/.9^2);
        h = h/sum(h);
        im = convn(im, h, 'same');
        im = convn(im, h', 'same');
        h = exp(-x.^2/2/1^2);
        h = h/sum(h);
        im = convn(im, reshape(h, 1,1,length(h)), 'same');
        for i=1:size(im,1)
            for j=1:size(im,2)
                if im(i,j)<stack(i,j,k)/2
                    im(i,j)=stack(i,j,k);
                end
                if stack(i,j,k)==0
                    im(i,j)=NaN;
                end
            end
        end
        
        
        % im(im==0)=NaN;
        % im=smooth2a(im,2,2);
        
        %gradient
        [fx, fy] = gradient(im*cline_para.gradient_force); %computing the gradient
        fx(isnan(fx))=0;
        fy(isnan(fy))=0;
        
        %% prepare for cline_para.iterations
        xys_avg = zeros(size(xys));
        counter = 0;
        %moving the snake in each iteration
        if k==1
            iterations=cline_para.iterations*6;
        else
            iterations=cline_para.iterations;
        end
        %%
        for i=1:iterations
            %      waitbar(((k-1)*cline_para.iterations+i+(k>1)*cline_para.iterations*9)/((9+length(time_idx))*cline_para.iterations));
            
            %calculate image gradient force on the nodes
            if cline_para.interpolate_gradient_flag
                % actual interpolation, slow!!!
                fs = [interp2(fx,xs,ys,'*linear'), interp2(fy,xs,ys,'*linear')];
                fs(isnan(fs)) = 0;
            else
                
                s_ind = round(xys)    ;
                fs = zeros(m, 2);
                %& ~isnan(im(sub2ind([61,61],s_ind(:,2),s_ind(:,1))))
                in_image_flag = all((s_ind>2)&(bsxfun(@le, s_ind, [col-2, row-2])), 2) ;
                gradient_ind = sub2ind([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
                in_image_flag(in_image_flag)=~isnan(im(gradient_ind)) & in_image_flag(in_image_flag);
                
                gradient_ind = sub2ind([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
                
                %
                % FX=fx(gradient_ind);
                % FY=fy(gradient_ind);
                %                 fs=-[FX,FY];
                fs=zeros(length(in_image_flag),2);
                fs(in_image_flag, :) = -[fx(gradient_ind), fy(gradient_ind)];
                
                
            end
            
            
            
            %viscos force between frames on all nodes:
            if k > 1 && pass ==1
                fs = fs - (xys_predicted(:,:,k-1)-xys) .* cline_para.inter_frame_viscos_factor .* 0.1;

            elseif k > 2 && pass == 1
                fs = fs + (xys_predicted(:,:,k-1)-xys) .* cline_para.inter_frame_viscos_factor .* 0.1; %.5

            elseif pass == 2 && k>1
                %    second round
                %   fs = fs + (xys_predicted(:,:,k)-xys) .* cline_para.inter_frame_viscos_factor .* 5;
                
            end
            
            %             size(fs)
            %             size(xys)
            %calculate the new position of snake
            f_xys=xys; %save old position
            new_xys = Ainv * (cline_para.gamma*xys - cline_para.kappa*fs);
            delta_xyso=new_xys-f_xys; %find change, ignore small changes
            %          delta_xyso(abs(delta_xyso)<.002)=0;
            
            if rem(i,cline_para.it_loop)<cline_para.group_iterations ||...
                    (k==1 && rem(i,cline_para.it_loop)<(cline_para.it_loop-10) && i<300)
                
                delta_xyso(~in_image_flag,:)=0;
                %               delta_c=min(max(delta_xyso))*mean(normr(delta_xyso));
                %      delta_c=min(max(abs(delta_xyso)))*mean(normr(delta_xyso));
                delta_c=mean(delta_xyso);
            else if k~=1 || i>cline_para.group_iterations
                    delta_c=mean(delta_xyso);
                else  if k==1 && i <cline_para.group_iterations
                        delta_c=[0,0];
                    end
                end
            end
            
            delta_xyso=[delta_xyso(:,1)-delta_c(1),delta_xyso(:,2)-delta_c(2)];
            if (k>1 && rem(i,cline_para.it_loop)> cline_para.group_iterations)||...
                    (k==1 && rem(i,cline_para.it_loop)> cline_para.group_iterations&& i>=300)
                
                delta_r=sum(delta_xyso.*n_vector,2);                %project steps radially and smooth around the circle
                h = exp(-x.^2/2/1^2);
                delta_r=repmat(delta_r,3,1);
                delta_r=conv(delta_r,h','same');
                delta_r=delta_r(cline_para.pts+1:2*cline_para.pts);
                
                %minradius check
                xys2=xys;
                xys2(~in_image_flag,:)=[];
                xys_d=sum(sqrt(sum(diff(xys2).^2,2)))/pi;
                if xys_d<minR %&& mean(delta_r)<-.01
                    
                    delta_r=(delta_r)*0;%&zeros(size(delta_r));
                    %             elseif xys_d >15 && mean(delta_r)>.01
                    %                 delta_r=-abs(delta_r);
                    %
                end
                
                
                
                
                delta_c=[0,0];
                %delta_r(all(fs==0,2))=0;
                
            else if rem(i,cline_para.it_loop)==cline_para.group_iterations
                    if k==1
                        delta_c=mean(xys)-[cline_para.center(1)+deltacount(1),cline_para.center(2)+deltacount(2)];
                    else if k>1
                            
                            xys2=xys;
                            xys2(~in_image_flag,:)=[];
                            xys_s=[0,cumsum(sqrt(sum(diff(xys2).^2,2)))'];
                            cx=mean(interp1(xys_s,xys2(:,1),1:max(xys_s)));
                            cy=mean(interp1(xys_s,xys2(:,2),1:max(xys_s)));
                            delta_c=[cx,cy]-[cline_para.center(1)+deltacount(1),cline_para.center(2)+deltacount(2)];
                        end
                    end
                end
                delta_r=zeros(length(n_vector),1);
                %  delta_c=mean(delta_xyso);
            end
            deltacount=deltacount+delta_c;
            
            
            delta_xys=zeros(size(delta_xyso));
            for ii=1:cline_para.pts
                delta_xys(ii,1:2)=n_vector(ii,1:2)*delta_r(ii)+delta_c;
                
            end
            
            
            %delta_xys=delta_xys-repmat(mean(n_vector(:,1:2).*repmat(delta_r,1,2)),length(delta_xys),1);
            %             delta_xys(:,1)=delta_xys(:,1)-mean(delta_xys(:,1));
            %             delta_xys(:,2)=delta_xys(:,2)-mean(delta_xys(:,2));
            
            xys=xys+delta_xys;
            
            if i > iterations -60 %smooth over the last 60 iterations
                
                counter = counter + 1;
                xys_avg = (xys_avg + xys);
            end
            
            
            
            
            
            if cline_para.show_flag==2
                
                if rem(i,20)==0
                    figure(1)
                    imagesc(im);
                    %caxis([min(im(im~=0)),max(max(im))])
                    axis equal
                    hold on;
                    quiver(xys(:,1),xys(:,2),fs(:,1),fs(:,2))
                    plot(xys(:,1), xys(:,2), 'r.-');
                    plot(xys_avg(:,1)/counter, xys_avg(:,2)/counter, 'b.-');
                    plot(xys_predicted(:,1, k), xys_predicted(:,2, k), 'c.-');
                    
                    %      scatter(mean(xys(:,1)), mean(xys(:,2)),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
                    scatter(cline_para.center(1)+deltacount(1),cline_para.center(2)+deltacount(2),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
                    hold off;
                    axis equal
                    pause(.005)
                    if cline_para.movie_flag==1;
                        frame = getframe;
                        writeVideo(writerObj,frame);
                    end
                    
                end
            end
            
            residue=(mean(xys));
            RR(:,:,i)=xys;
            if i>80
                %             RRR(i,:)=residue-RR(i-80,:);
                %             [i,mean(mean(RRR(i-60:i)))]
                RRR(:,:,i)=xys-RR(:,:,i-80);
                XYS_avg(:,:,i)=xys;
                if i>160
                    %                                 figure(2)
                    %                 plot(RRR(80:end,:))
                    if (k>1 && all(all(all(abs(RRR(:,:,i-40:i))<cline_para.threshold)))) ||i==iterations...
                            ||(k==1 && all(all(all(abs(RRR(:,:,i-40:i))<cline_para.threshold/10))))
                        xys_avg=sum(XYS_avg(:,:,i-40+1:i),3)/40;
                        xys_all(:,:,k) = xys_avg;
                        clear RRR
                        
                        break
                    end
                end
            end
            
            
            
            
            
        end
        
        if k == 1 && pass == 1
            xys_predicted(:,:,1) = xys_avg;
            xys_predicted(:,:,2) = xys_avg;
        elseif k < length(time_idx) && pass == 1
            xys_predicted(:,:,k+1) = xys_avg;
        end
        
        if k==1
            c1(2)=xys_avg(cline_para.pts,2);
            c1(1)=xys_avg(cline_para.pts*.75,1);
        end
        if k==stack_z_size
            c2(2)=xys_avg(cline_para.pts,2);
            c2(1)=xys_avg(cline_para.pts*.75,1);
        end
        if k == 1 && pass == 1
            xys_predicted(:,:,1) = xys_avg;
            xys_predicted(:,:,2) = xys_avg;
        elseif k < length(time_idx) && pass == 1
            xys_predicted(:,:,k+1) = xys_avg;
        end
        
        %         if pass == 2
        %             figure(2)
        %             cla;
        %          imshow(stack_z_size_max(:,:,k),[]);
        %          hold on;
        %          plot(xyzs(:,1), xyzs(:,2), 'r.-');
        %          hold off;
        %             writeVideo(vidObj,getframe);
        %         end
        % subplot(1,2,1)
        
        if cline_para.show_flag==1
            imagesc(im);
            axis equal
            hold on
            plot(xys_avg(:,1), xys_avg(:,2), 'g.-');
            %  plot(xys_predicted(:,1, k), xys_predicted(:,2, k), 'c.-');
            scatter(cline_para.center(1)+deltacount(1),cline_para.center(2)+deltacount(2),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
            caxis([min(min(stack(stack~=0))),max(max(stack(:,:,k)))])
            hold off
        end
        
        %   caxis([min(min(stack(stack~=0))),max(max(stack(:,:,k)))])
        
        
        
        
        
        %plot(xys_predicted(:,1,k), xys_predicted(:,2,k), 'b.-');
        %      scatter(mean(xys(:,1)),mean(xys(:,2)),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
        
        if cline_para.movie_flag==1 && pass==2;
            frame = getframe;
            writeVideo(writerObj,frame);
        end
        %    pause(.1);
        
        
        
    end
    % smooth predicted position for the 2nd pass
    xys_all_pad = padarray(xys_all, [0 0 2], 'replicate');
    xys_predicted = convn(xys_all_pad, reshape([1/5, 1/5, 1/5, 1/5, 1/5], 1,1,5), 'same');
    xys_predicted = xys_predicted(:,:,3:end-2);
    %xys_predicted=xys_all;
    % if k == length(time_idx);
    %     xys=xys_predicted(:,:,1);
    % end
    
    if all(all(fs==0)) ||abs(xys_avg(1,1))>2*row
        break
    end
    
    
end
% close(vidObj);
% save data:
%save('xys_all.mat', 'xys_all');
clearvars -except xys_all hyper_stack stack_z_size_max stack_t_size c1 c2 RR;