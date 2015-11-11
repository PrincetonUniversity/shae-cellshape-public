function xys_all = f_CircDet3(stack, cline_para,initial)

        cline_para.group_iterations=00;
        cline_para.iterations=100;
        
initial=squeeze(initial);
%stack=permute(stack,[3,1,2]);
iterations=cline_para.iterations;
deltacount=[0,0];
    
if size(size(stack))==3
    stack1=zeros(size(stack),1);
    stack1(:,:,1)=stack;
    stack=stack1;
end
[row, col, stack_z_size] = size(stack);
kernal=fspecial('sobel');





time_idx = 1:stack_z_size;

% initialize snake

theta_interp = -pi+2*pi/cline_para.pts:2*pi/cline_para.pts:pi;
n_vector=zeros(cline_para.pts,2);
n_vector(:,1)=cos(theta_interp);
n_vector(:,2)=sin(theta_interp);



xys=initial;
center=mean(xys);


anchor1=xys(round(cline_para.pts/2),:);
anchor2=xys(cline_para.pts,:);


[m n] = size(xys);
    
% populating the penta diagonal matrix
A = zeros(m,m);
brow = zeros(1,m);
brow(1,1:5) = [cline_para.beta, -(cline_para.alpha + 4*cline_para.beta), (2*cline_para.alpha + 6 *cline_para.beta), -(cline_para.alpha + 4*cline_para.beta), cline_para.beta];
brow=circshift(brow',-2)';
% A(1, 1:3) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];
% A(2, 1:4) = [-cline_para.alpha-2*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-4*cline_para.beta, cline_para.beta];
for i=1:m
    A(i,:) = brow;
    brow = circshift(brow',1)'; % Template row being rotated to egenrate different rows in pentadiagonal matrix
end
% A(m-1, m-3:m) = [cline_para.beta, -cline_para.alpha-4*cline_para.beta, 2*cline_para.alpha+5*cline_para.beta, -cline_para.alpha-2*cline_para.beta];
% A(m, m-2:m) = [cline_para.beta, -2*cline_para.beta, cline_para.beta];


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
for pass = 1
    %k: stack_z_size index
    for k = 1:length(time_idx)
%         if k==1
%             cline_para.iterations=1;
%         else
%             cline_para.iterations=iterations;
%         end
        
        goodflag=0;
        %obtain stack_z_size(t)
        im = stack(:,:,k); 
        %% filter image to "normalize" intensity

%         if cline_para.gradflag==1;
%             imMask=(im~=0);
%             imMask=smooth2a(imMask,2,2);
%             imMask=(imMask==1);
%             hoIm=conv2(im,kernal,'same');
%             vertIm=conv2(im,kernal','same');
%             im=sqrt(hoIm.^2+vertIm.^2);
%             im=imMask.*im;
%             
%                     im(1:2,:)=0;
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
            
%             % heavy blur
%             h = exp(-x.^2/2/3.5^2);  
%             h = h/sum(h);
%             im = convn(im, h, 'same');
%             im = convn(im, h', 'same');
%             h = exp(-x.^2/2/2.5^2);
%             h = h/sum(h);
%             im = convn(im, reshape(h, 1,1,length(h)), 'same');
%            % filter
%             im = im1 ./im; 
%         else
%             x = -30:30;
%             h = exp(-x.^2/2/10^2); 
%             h = h/sum(h);
%             im = convn(im, h, 'same');
%             im = convn(im, h', 'same');
%             im = -im/3000;
%         end
        
        
        
        

        
        
        
                    % light blur
                    
                    
                    
%             x = -12:12;
%             h = exp(-x.^2/2/1.5^2); 
%             h = h/sum(h);
%             im = convn(im, h, 'same');
%             im = convn(im, h', 'same');
%             h = exp(-x.^2/2/1^2); 
%             h = h/sum(h);
%             im1 = convn(im, reshape(h, 1,1,length(h)), 'same');
%             % heavy blur
%             h = exp(-x.^2/2/3.5^2);  
%             h = h/sum(h);
%             im = convn(im, h, 'same');
%             im = convn(im, h', 'same');
%             h = exp(-x.^2/2/2.5^2);
%             h = h/sum(h);
%             im = convn(im, reshape(h, 1,1,length(h)), 'same');
%             %filter
%             im = im1 ./im;
%         
%         
%         
        
        
        
        %gradient
        [fx, fy] = gradient(im*cline_para.gradient_force); %computing the gradient
        fx(isnan(fx))=0;
fy(isnan(fy))=0;
        %% prepare for cline_para.iterations
        xys_avg = zeros(size(xys));
        counter = 0;
        %moving the snake in each iteration
        for i=1:cline_para.iterations
            if goodflag==0;
            %calculate image gradient force on the nodes
            if cline_para.interpolate_gradient_flag
                % actual interpolation, slow!!!
                fs = [interp2(fx,xs,ys), interp2(fy,xs,ys)];
                fs(isnan(fs)) = 0;
            else
                s_ind = round(xys);        
                fs = zeros(m, 2);    
                in_image_flag = all((s_ind>2)&(bsxfun(@le, s_ind, [col-2, row-2])), 2);
                gradient_ind = sub2ind([row, col], s_ind(in_image_flag, 2), s_ind(in_image_flag, 1));
                FX=interp2(fx,xys(:,1),xys(:,2),'*linear',0);
                FY=interp2(fy,xys(:,1),xys(:,2),'*linear',0);
                
                %fs(in_image_flag, :) = -[fx(gradient_ind), fy(gradient_ind)];
                fs=-[FX,FY];
            end
           
            fs(cline_para.pts,:)=[0,0];
            fs(round(cline_para.pts/2),:)=[0,0];
            
           
%           viscos force between frames on all nodes:        
%             if k > 1 && pass ==1
%                 fs = fs + (xys_predicted(:,:,k)-xys) .* cline_para.inter_frame_viscos_factor .* 0.1;
%             elseif k > 2 && pass == 1
%                 fs = fs + (xys_predicted(:,:,k)-xys) .* cline_para.inter_frame_viscos_factor .* 0.5;
%             elseif pass == 2
%                 %second round
%                 fs = fs + (xys_predicted(:,:,k)-xys) .* cline_para.inter_frame_viscos_factor .* 5;
%             end

            
            %calculate the new position of snake
            f_xys=xys;
            new_xys = Ainv * (cline_para.gamma*xys - cline_para.kappa*fs);
            delta_xyso=new_xys-f_xys;
            
            if k~=1 && i<cline_para.group_iterations

                delta_c=min(max(delta_xyso))*mean(normr(delta_xyso(1:round(cline_para.pts/5):cline_para.pts,:)));
                
            else if k~=1 || i>cline_para.group_iterations
                    delta_c=mean(delta_xyso);
                else  if k==1
                        delta_c=[0,0];
                    end
                end
            end
            
            
            
            
            
            
            delta_xyso=[delta_xyso(:,1)-delta_c(1),delta_xyso(:,2)-delta_c(2)];
    if i>    cline_para.group_iterations||k==1
            delta_r=sum(delta_xyso.*n_vector,2);
            delta_c=[0,0];
    else if i==cline_para.group_iterations
            delta_c=mean(xys)-[center(1)+deltacount(1),center(2)+deltacount(2)];
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

            change_xys=f_xys/cline_para.gamma-xys;
            % take average value at the end of the cline_para.iterations
           changes=(sqrt(change_xys(:,1).^2+change_xys(:,2).^2));
           
 
                                       xys(cline_para.pts,:)=anchor2;
           xys(round(cline_para.pts/2),:)=anchor1;           
            if i > cline_para.iterations - 60
               
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
                quiver(xys(:,1),xys(:,2),delta_xyso(:,1),delta_xyso(:,2))
                 plot(xys(:,1), xys(:,2), 'r.-');
                 plot(initial(:,1),initial(:,2),'g.-');
                %plot(xys_avg(:,1)/counter, xys_avg(:,2)/counter, 'b.-');
          %      scatter(mean(xys(:,1)), mean(xys(:,2)),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)     
                scatter(cline_para.center(1)+deltacount(1),cline_para.center(2)+deltacount(2),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)                
                hold off;
                axis equal
                pause(.005)

                end
            end




             end
            
        end
        xys_avg = xys_avg/ counter ;
        xys_all(:,:,k) = xys_avg;
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
       %imagesc(stack(:,:,k)); 
       imagesc(im)
              axis equal
        hold on
   %     plot(xys_predicted(:,1, k), xys_predicted(:,2, k), 'c.-');
        plot(xys_avg(:,1), xys_avg(:,2), 'g.-');
        scatter(center(1)+deltacount(1),center(2)+deltacount(2),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
        scatter(mean(xys(:,1)),mean(xys(:,2)),'DisplayName','xys_all(:,2) vs. xys_all(:,1)','YDataSource','xys_all(:,2)');figure(gcf)
        hold off
end
  %    pause(.1);
    end
    % smooth predicted position for the 2nd pass
    xys_all_pad = padarray(xys_all, [0 0 2], 'replicate');
    xys_predicted = convn(xys_all_pad, reshape([1/5, 1/5, 1/5, 1/5, 1/5], 1,1,5), 'same');
    xys_predicted = xys_predicted(:,:,3:end-2);
end
% close(vidObj);

% save data:
%save('xys_all.mat', 'xys_all');
clearvars -except xys_all hyper_stack stack_z_size_max stack_t_size;