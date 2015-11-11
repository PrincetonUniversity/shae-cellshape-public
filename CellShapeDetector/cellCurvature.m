function [xf,yf,zf,gc,gm,kplus,kminus]=cellCurvature(x,y,z,Nv,Bv)
%this program takes cell coordinate matricies x,y and z, along with
%centerline Frenet vectors of the centerline to make a curvature map over
%all points. The original coordinate system of the cell is one that is
%cylindrical for the body and polar at the poles, resulting in a
%singularity at each pole. For calculating curvatures at the poles, the
%first and last end_slices are transfered into polar projection coordinates
%using the midpoint "ender"slices from the pole.


if size(x,2)>15
    end_slices=14;
    ender=15;
else
    ender=size(x,2);
    end_slices=size(x,2)-1;
end

start_pt=mean([x(:,ender),y(:,ender),z(:,ender)]);
end_pt=mean([x(:,end-ender+1),y(:,end-ender+1),z(:,end-ender+1)]);
xps=x(:,1:end_slices)-start_pt(1);
yps=y(:,1:end_slices)-start_pt(2);
zps=z(:,1:end_slices)-start_pt(3);

xpe=x(:,end-end_slices+1:end)-end_pt(1);
ype=y(:,end-end_slices+1:end)-end_pt(2);
zpe=z(:,end-end_slices+1:end)-end_pt(3);

xps=cat(1,xps,xps(1,:));
yps=cat(1,yps,yps(1,:));
zps=cat(1,zps,zps(1,:));
xpe=cat(1,xpe,xpe(1,:));
ype=cat(1,ype,ype(1,:));
zpe=cat(1,zpe,zpe(1,:));
[X,Y]=meshgrid(1:.5:size(xps,2),(1:.5:size(xps,1))');
xps=interp2(xps,X,Y,'*spline');
yps=interp2(yps,X,Y,'*spline');
zps=interp2(zps,X,Y,'*spline');
xpe=interp2(xpe,X,Y,'*spline');
ype=interp2(ype,X,Y,'*spline');
zpe=interp2(zpe,X,Y,'*spline');

poleSize=numel(xps);
xpsl=reshape(xps,poleSize,1);
ypsl=reshape(yps,poleSize,1);
zpsl=reshape(zps,poleSize,1);
xpel=reshape(xpe,poleSize,1);
ypel=reshape(ype,poleSize,1);
zpel=reshape(zpe,poleSize,1);

Ns=xpsl*Nv(1,1)+ypsl*Nv(1,2)+zpsl*Nv(1,3);
Ne=xpel*Nv(end,1)+ypel*Nv(end,2)+zpel*Nv(end,3);
Bs=xpsl*Bv(1,1)+ypsl*Bv(1,2)+zpsl*Bv(1,3);
Be=xpel*Bv(end,1)+ypel*Bv(end,2)+zpel*Bv(end,3);

warning('off','all');

Fxs=TriScatteredInterp(Ns,Bs,xpsl);
Fxe=TriScatteredInterp(Ne,Be,xpel);
Fys=TriScatteredInterp(Ns,Bs,ypsl);
Fye=TriScatteredInterp(Ne,Be,ypel);
Fzs=TriScatteredInterp(Ns,Bs,zpsl);
Fze=TriScatteredInterp(Ne,Be,zpel);

warning('on','all');

[qx,qy]=meshgrid(-1000:20:1000,-1000:20:1000);

xs=Fxs(qx,qy);
ys=Fys(qx,qy);
xe=Fxe(qx,qy);
ye=Fye(qx,qy);
zs=Fzs(qx,qy);
ze=Fze(qx,qy);
[gme,gce]=MGCurve(xe,ye,ze,5);
[gms,gcs]=MGCurve(xs,ys,zs,5);
gms=sign(mean(gms(~isnan(gms)))).*gms;
gme=sign(mean(gme(~isnan(gme)))).*gme;
gms=medfilt2(gms,[5,1]);
gme=medfilt2(gme,[5,1]);




gms=smooth2a(gms,10,10);
gcs=smooth2a(gcs,10,10);
gme=smooth2a(gme,10,10);
gce=smooth2a(gce,10,10);

Ns=reshape(Ns,size(xps));
Bs=reshape(Bs,size(xps));
Ne=reshape(Ne,size(xps));
Be=reshape(Be,size(xps));

gmsf=interp2(qx,qy,gms,Ns(1:2:end,1:2:end),Bs(1:2:end,1:2:end),'*linear');
gmef=interp2(qx,qy,gme,Ne(1:2:end,1:2:end),Be(1:2:end,1:2:end),'*linear');
gcef=interp2(qx,qy,gce,Ne(1:2:end,1:2:end),Be(1:2:end,1:2:end),'*linear');
gcsf=interp2(qx,qy,gcs,Ns(1:2:end,1:2:end),Bs(1:2:end,1:2:end),'*linear');




[gm,gc,xf,yf,zf]=MGCurveloop(x,y,z,5);
%replace gc and gm at ends with curvatures from original
%paramaterization

gm=sign(median(gm(~isnan(gm)))).*gm;



%         end_change=min(floor(size(x,2)/3),floor(size(gcsf,2)*2/3));
%     %    end_change=min(end_change,8);
end_change=max(sum(~isnan(sum(gcsf))),sum(~isnan(sum(gcef))));
for i=1:end_change
    if nanstd(gc(:,i))>nanstd(gcsf(:,i)) || isnan(sum(gc(:,i)))
        gc(:,i)=gcsf(:,i);
        gm(:,i)=gmsf(:,i);
    end
    if nanstd(gc(:,end+1-i))>nanstd(gcef(:,end+1-i)) || isnan(sum(gc(:,end-i+1)))
        gc(:,end-i+1)=gcef(:,end-i+1);
        gm(:,end-i+1)=gmef(:,end-i+1);
    end
end



%calculate other curvatures, -0 coordinates mean only body points
kplus=gm+sqrt(max(gm.^2-gc,0));
kminus=gm-sqrt(max(gm.^2-gc,0));

