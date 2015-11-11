 function [P_start,P_end,theta_start,theta_end, Bp_start,Bp_end,start_pt,end_pt]=cellPoleCoordSys(hyper_stack,Nv,Bv,c1,c2,xyzs,pixel_interp,window,flag,end_angles)
    % this function creates a stack of the poles of the cell parameterized
    % by rotating the Bv vectors around the Nv vectors at the poles about
    % the points c1 and c2 and taking slices.
 %number of slice angles to fit for poles
    xpole_slice=zeros(window*pixel_interp+1,2*window*pixel_interp+1,end_angles);
    ypole_slice=xpole_slice;
    zpole_slice=xpole_slice;
    Bp_start=zeros(end_angles,3);
   Bp_end=zeros(end_angles,3);
    theta_start=-(0:end_angles-1)*pi/(end_angles-1);
    theta_end=-theta_start;
    Bp=zeros(end_angles,3);
    
    %create basis for poles at start and end by rotating B vectors around N
    %by theta
    for i=1:1:end_angles
        Bp_start(i,:)=rodrigues_rot(Bv(1,:),Nv(1,:),theta_start(i));
        Bp_end(i,:)=rodrigues_rot(Bv(end,:),Nv(end,:),theta_end(i));
    end
    %go through 2 runs, one for start and one for end
    for run=1:2 %set up for each pole
        if flag.endpole==1 && run==1
            Np=Nv(end,:);
            Bp1=Bv(end,:);
            theta=theta_end;
            c2f=c2*[Np',Bp1']';
            end_pt=xyzs(end,:)+c2f;

            go=1;
        elseif flag.endpole==0 && run==1
                go=0;
        end


        if flag.startpole==1 && run==2
            Np=Nv(1,:);
            Bp1=Bv(1,:);
            theta=theta_start;
            c1f=c1*[Np',Bp1']';
            start_pt=xyzs(1,:)+c1f;

        elseif flag.startpole==0 && run==2
                go=0;
            
        end
        
        if go==1;
            
            [J,K]=meshgrid(-window:1/pixel_interp:0,...
                -window:1/pixel_interp:window);
            xpole_slice1=J*Bp1(1)+K*Np(1);
            ypole_slice1=J*Bp1(2)+K*Np(2);
            zpole_slice1=J*Bp1(3)+K*Np(3);
            
            
            for i=1:1:end_angles
                Bp(i,:)=rodrigues_rot(Bp1,Np,theta(i));
                
                for k=1:pixel_interp*2*window+1
                    temp=rodrigues_rot([xpole_slice1(k,:)',ypole_slice1(k,:)',zpole_slice1(k,:)'],Np,theta(i));
                    
                    if run==1
                        xpole_slice(:,k,i)=temp(:,1)+end_pt(1);
                        ypole_slice(:,k,i)=temp(:,2)+end_pt(2);
                        zpole_slice(:,k,i)=temp(:,3)+end_pt(3);
                    elseif run==2
                        xpole_slice(:,k,i)=temp(:,1)+start_pt(1);
                        ypole_slice(:,k,i)=temp(:,2)+start_pt(2);
                        zpole_slice(:,k,i)=temp(:,3)+start_pt(3);
                        
                    end
                    
                end
            end
        end
        
        %interpolate to find intensities in this basis
        
%         if flag.gradient==0
%             P=interp3(hyper_stack,ypole_slice,xpole_slice,zpole_slice,'*linear',0);
%         elseif flag.gradient==1
%             P=interp3(grad_stack,ypole_slice,xpole_slice,zpole_slice,'*linear',0);
%             
%         end
        
                            P=interp3(hyper_stack,ypole_slice,xpole_slice,zpole_slice,'*linear',0);

        if flag.endpole==1 && run==1
            P_end=[P(1:pixel_interp*window,:,:)
                flipdim(P,1)];
        end
        
        
        if flag.startpole==1 && run==2
            P_start=[P(1:pixel_interp*window,:,:)
                flipdim(P,1)];
        end
        
    end
    