function [xp,yp,zp,tip]=cellPoleCoordSys2Image(xys_end,Np,Bp,cline_para,end_angles,pixel_interp)
%takes contours xys_end in coordinate system Np,Bp and transforms back into
%cartesian coordinates.
           % Np=Nv(plane_num,:);

            Tv=cross(Np,Bp(1,:));
            xyz_end_final=zeros(round(cline_para.pts/2)+1,3,end_angles);
            
            for i=1:end_angles
                xyz_end= [xys_end(end,:,i)',xys_end(1:round(cline_para.pts/2),:,i)']'*[Np',Bp(i,:)']'/pixel_interp;
                xyz_end_final(:,:,i)=xyz_end;
            end
            xyz_end_final(:,1:2,:)=xyz_end_final(:,1:2,:);
            xyz_end_final(:,3,:)=xyz_end_final(:,3,:);
            
            xp=squeeze(xyz_end_final(:,1,:))';
            yp=squeeze(xyz_end_final(:,2,:))';
            zp=squeeze(xyz_end_final(:,3,:))';
            tip=max(max(abs(xp*Tv(1,1)+yp*Tv(1,2)+zp*Tv(1,3))));
