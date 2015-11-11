function [xps,yps,zps,end_slices]=PoleRotate(l_end,cline_para,pv,xp,yp,zp,theta)

% Fix Ends and reparameterize into coordinate system that matches the body
    %This converts pole points into spherical coordinates, rotates the
    %coordinate axis to a new phi and theta, and then fixes the angle to
    %replace it at the poles. It also calculates curvatures in the original
    %coordinate system to use later for replacing the singularity at the
    %the tips of the cell.
    
    end_slices=ceil(l_end*pi/2);
    alphaE=pi/2:-pi/(2*end_slices):0;
    beta=-pi/2:2*pi/cline_para.pts:3*pi/2;
    
  %create angle basis
    [BETAE,ALPHAE]=meshgrid(beta,alphaE);
    %this is the transformation for the rotation of the spherical
    %coordinate system, alpha and theta2 are the old and new polar
    %angles
    
    THETA2E=-acos(sin(ALPHAE).*sin(BETAE));
    CAE=cot(ALPHAE);
    CAE(abs(cos(BETAE(:,1)))<.001,1)=0;
    PHI2E=atan(-CAE./cos(BETAE));

    
    
    phi=-pi:2*pi/cline_para.pts:0;
    phi=-phi;
    
        
        
            PHI2=atan(-CAE./cos(BETAE));
            PHI2(PHI2<0)=pi+PHI2(PHI2<0);
            PHI2(:,1)=pi/2;
            PHI2=pi-PHI2;
            xyzp=[xp(1,1),yp(1,1),zp(1,1)];
            
            %PHI2=flipdim(PHI2,2);
            
            
            [THETA,PHI]=meshgrid(phi,theta);
           
                   FX=scatteredInterpolant(THETA(:),PHI(:),xp(:),'natural','none');
                   FY=scatteredInterpolant(THETA(:),PHI(:),yp(:),'natural','none');
                   FZ=scatteredInterpolant(THETA(:),PHI(:),zp(:),'natural','none');

                   PHI2E=mod(PHI2E,pi);
if all(theta<=0)
PHI2E=-PHI2E;
end
if all(phi>=0);
    THETA2E=-THETA2E;
end

xps=FX(THETA2E,PHI2E);
yps=FY(THETA2E,PHI2E);
zps=FZ(THETA2E,PHI2E); 
   

        
            xps=xps(:,2:end);
            yps=yps(:,2:end);
            zps=zps(:,2:end);
           