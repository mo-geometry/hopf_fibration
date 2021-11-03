function plot2sphere(hopf,nCircles)

% initialize
[~,~,nFrames] = size(hopf.Q);
[ sph.x, sph.y, sph.z ] = sphere();
% generate sphere figure

% plot points
for ii=1:nFrames
    
    
    
    %% plot
    figure(1);clf
    subplot(3,2,[5,6])
    axis equal;     grid off;    hold on;   axis off
    surf( sph.x,sph.y,sph.z); 
    alpha(0.5)
    view([30 30])
    colormap gray
    shading interp
    for jj=1:nCircles
        if isnan(hopf.pts(jj,1,ii)), continue; end
%         colorV = [hopf.pts(jj,1,ii)*0.5+0.5,hopf.pts(jj,2,ii).^2,hopf.pts(jj,3,ii).^2];
        colorV = [hopf.pts(jj,1,ii)*0.5+0.5,hopf.pts(jj,2,ii)*0.5+0.5,hopf.pts(jj,3,ii)*0.5+0.5];
        plot3(hopf.pts(jj,3,ii),hopf.pts(jj,2,ii),hopf.pts(jj,1,ii),'*','color',colorV,'LineWidth',5)     
    end
    subplot(3,2,[1,2,3,4])
    axis equal;     grid off;    hold on;   axis off
    axis([-3 3 -3 3 -3 3])
    view([-30 15])
    for jj=1:nCircles
        if isnan(hopf.pts(jj,1,ii)), continue; end
%         disp(hopf.fibers(:,:,jj,ii))
        plot3(hopf.fibers(:,3,jj,ii),hopf.fibers(:,2,jj,ii),hopf.fibers(:,1,jj,ii),'-','color',hopf.pts(jj,:,ii).^2,'LineWidth',2)
%         plot3(hopf.fibers(:,3,jj,ii),hopf.fibers(:,2,jj,ii),hopf.fibers(:,1,jj,ii),'r-','LineWidth',2)
    end
    
end
