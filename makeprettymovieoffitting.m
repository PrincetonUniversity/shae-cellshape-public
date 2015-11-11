function makeprettymovieoffitting(savedSurface);
%% assumptions
PIXEL_SIZE = 0.08;% pixel size in microns
PIX_OFFSET = 0.3/PIXEL_SIZE; % offset for pseudo origin RGB vector base
%% find the maximal window size
% extract all of the surface points
allFs = [savedSurface.F];
allPoints = cell2mat({allFs.vertices}');

% find their range
xRange(1) = min(allPoints(:,1))-PIX_OFFSET;
xRange(2) = max(allPoints(:,1))+PIX_OFFSET;

yRange(1) = min(allPoints(:,2))-PIX_OFFSET;
yRange(2) = max(allPoints(:,2))+PIX_OFFSET;

zRange(1) = min(allPoints(:,3))-PIX_OFFSET;
zRange(2) = max(allPoints(:,3))+PIX_OFFSET;

%% find the colorscale for forcing

allForces = cell2mat({savedSurface.fs}');
colorrange(1) = -1* max(sqrt(sum(allForces.^2,2)));
colorrange(2) = +1*max(sqrt(sum(allForces.^2,2)));

%% loop over each image
[xx,yy,zz] = ndgrid(xRange,yRange,zRange);
figure(gcf);
set(gcf,'color','w');


for iiStep = 1
    subplot(3,1,1);
    cla
    scatter3(xx(:),yy(:),zz(:),'wx');
    hold on;
    
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET/2,...
        0.5/PIXEL_SIZE,0,0,...
        0,'r','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET/2,...
        0,0.5/PIXEL_SIZE,0,...
        0,'g','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET/2,...
        0,0,0.5/PIXEL_SIZE,...
        0,'b','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    
    p1 = patch(allFs(iiStep),'facevertexCdata',sqrt(sum((savedSurface(iiStep).fs.^2),2)),'facecolor','interp','facealpha',0.8,'edgecolor','k','edgealpha',0.15);
    set(gca,'Clim',colorrange);
    xlim(xRange);
    ylim(yRange);
    zlim(zRange);
    set(gca,'DataAspectRatio',[1,1,1]);
    view(0,90);
    figure(gcf);
    axis off
    
    lighting gouraud
    camlight('right')
    material dull
    pause(0.1);
    
    subplot(3,1,2);
    cla
    scatter3(xx(:),yy(:),zz(:),'wx');
    hold on;
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0.5/PIXEL_SIZE,0,0,...
        0,'r','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0,0.5/PIXEL_SIZE,0,...
        0,'g','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0,0,0.5/PIXEL_SIZE,...
        0,'b','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    
    
    p2 = patch(allFs(iiStep),'facevertexCdata',sqrt(sum((savedSurface(iiStep).fs.^2),2)),'facecolor','interp','facealpha',0.8,'edgecolor','k','edgealpha',0.15);
    set(gca,'Clim',colorrange);
    xlim(xRange);
    ylim(yRange);
    zlim(zRange);
    set(gca,'DataAspectRatio',[1,1,1]);
    view(0,0);
    figure(gcf);
    axis off
    lighting gouraud
    camlight('right')
    material dull
    pause(0.1);
    
    subplot(3,1,3);
    cla
    scatter3(xx(:),yy(:),zz(:),'wx');
    hold on;
    
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0.5/PIXEL_SIZE,0,0,...
        0,'r','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0,0.5/PIXEL_SIZE,0,...
        0,'g','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    h =  quiver3(xRange(1)+PIX_OFFSET/2,yRange(1)+PIX_OFFSET/2,zRange(1)+PIX_OFFSET,...
        0,0,0.5/PIXEL_SIZE,...
        0,'b','linewidth',3);
    adjust_quiver_arrowhead_size(h,8);
    
    p3 = patch(allFs(iiStep),'facevertexCdata',sqrt(sum((savedSurface(iiStep).fs.^2),2)),'facecolor','interp','facealpha',0.8,'edgecolor','k','edgealpha',0.15);
    set(gca,'Clim',colorrange);
    xlim(xRange);
    ylim(yRange);
    zlim(zRange);
    set(gca,'DataAspectRatio',[1,1,1]);
    view(90,0);
    figure(gcf);
    axis off
    lighting gouraud
    camlight('right')
    material dull
    pause(0.1);
    
end

for iiStep = 2:length(savedSurface);
    subplot(3,1,1);
    delete(p1);
    
    p=allFs(iiStep).vertices;
    tri=allFs(iiStep).faces;
    normvec=vertnorm(tri,p);
    
    p1 = patch(allFs(iiStep),'facevertexCdata',-1*dot(normvec,savedSurface(iiStep).fs,2),'facecolor','interp','facealpha',0.85,'edgecolor','k','edgealpha',0.03);
    material dull
    %     set(p1,'AmbientStrength',0.6,'DiffuseStrength',0.7,'SpecularStrength',0.3);
    
    subplot(3,1,2);
    delete(p2)
    p2 = patch(allFs(iiStep),'facevertexCdata',-1*dot(normvec,savedSurface(iiStep).fs,2),'facecolor','interp','facealpha',0.85,'edgecolor','k','edgealpha',0.03);
    material dull
    %     set(p2,'AmbientStrength',0.6,'DiffuseStrength',0.7,'SpecularStrength',0.3);
    
    subplot(3,1,3);
    delete(p3)
    p3 = patch(allFs(iiStep),'facevertexCdata',-1*dot(normvec,savedSurface(iiStep).fs,2),'facecolor','interp','facealpha',0.85,'edgecolor','k','edgealpha',0.03);
    material dull
    %     set(p3,'AmbientStrength',0.6,'DiffuseStrength',0.7,'SpecularStrength',0.3);
    export_fig(['saveingAMovie/plane',num2str(iiStep)], '-tif', '-nocrop',gcf);
    
end



