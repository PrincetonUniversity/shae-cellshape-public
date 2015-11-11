function xyzs_avg = ActiveContourFit2(im, cline_para, cline_initial)

if isempty(cline_para);
    cline_para.iterations=700;
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
cline_para.gamma=10; %step size
cline_para.gradflag=0;
cline_para.gradient_force =.5;%5
cline_para.interpolate_gradient_flag =0;
cline_para.show_flag=2; %1 to show end of every frame, 2 to show every compelted slice
cline_para.movie_flag=0;
end

if isempty(cline_initial)
    cline_initial=[repmat((size(im,2)/2),size(im,1),1),(1:size(im,1))'];
end






[row col, stack_z_size] = size(im);

% initialize snake
xyzs = cline_initial;
%s_interp = 0:0.1:s(end);
%xyzs = interp1(s, cline_initial, s_interp, 'spline');


[m n] = size(xyzs);
    
% populating the penta diagonal matrix
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
[L U] = lu(A + cline_para.gamma .* eye(m,m));
Ainv = inv(U) * inv(L); % Computing Ainv using LU factorization




        %% gradient
        
        [fy, fx, fz] = gradient(im*cline_para.gradient_force); %computing the gradient
        %% prepare for cline_para.iterations
        xyzs_avg = zeros(size(xyzs));
        counter = 0;
        %moving the snake in each iteration
        for i=1:cline_para.iterations;
            %calculate image gradient force on the nodes
             % actual interpolation, slow!!!
                xs=xyzs(:,1);
                ys=xyzs(:,2);
                zs=xyzs(:,3);
                fs = [interp3(fx,ys,xs, zs,'*linear'),...
                    interp3(fy,ys,xs, zs,'*linear'),...
                    interp3(fz,ys,xs, zs,'*linear')];
                
                fs(isnan(fs)) = 0;

            
            %confine the snake in the image
    %        fs(:, 3) = fs(:, 3) + (1./(1+exp((xyzs(:,3)+0)/3))-1./(1+exp((10-xyzs(:,3))/3)))*cline_para.z_confine_factor;
            % viscos force between frames on all nodes:        
xyzstart=xyzs(1,:);
xyzend=xyzs(end,:);
            % stretch or shrink the snake at the ends            
%             if cline_para.stretch_ends_flag
%                 s_head = xyzs(1,:) - xyzs(4, :);
%                 s_head = s_head/norm(s_head);
%                 fs(1, :) = fs(1, :) + s_head .* cline_para.stretching_force_factor(1);
%                 s_tail = xyzs(end, :) - xyzs(end - 3, :);
%                 s_tail = s_tail/norm(s_tail);
%                 fs(end, :) = fs(end, :) + s_tail .* cline_para.stretching_force_factor(2);
%             end
            %calculate the new position of snake
            xyzs = Ainv * (cline_para.gamma*xyzs + cline_para.kappa*fs);
            % take average value at the end of the cline_para.iterations
            if i > cline_para.iterations - 50
                xyzs_avg = xyzs_avg + xyzs;
                counter = counter + 1;
            end
            
            %end points are fixed
            xyzs(1,:)=xyzstart;
            xyzs(end,:)=xyzend;
            L=length(xyzs)-1;
            %redistribute points evenly
s=[0,cumsum(sqrt(diff(xyzs(:,1)).^2+diff(xyzs(:,2)).^2+diff(xyzs(:,3)).^2))']';
 
    xyzs(:,1)=interp1(s,xyzs(:,1),0:max(s)/L:max(s),'spline');
    xyzs(:,2)=interp1(s,xyzs(:,2),0:max(s)/L:max(s),'spline');
    xyzs(:,3)=interp1(s,xyzs(:,3),0:max(s)/L:max(s),'spline');
          
    %view contour
    
    scatter3(xyzs(:,1),xyzs(:,2),xyzs(:,3));
          axis equal
    
          pause(.1)
          
        end
         %avg last 50 iterations
         xyzs_avg = xyzs_avg / counter;

end
