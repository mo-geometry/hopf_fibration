function plot2sphere3(hopfA,hopfB,hopfC,nCircles)

% initialize
[~,~,nFrames] = size(hopfA.Q);
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
        if isnan(hopfA.pts(jj,1,ii)), continue; end
%         colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii).^2,hopfA.pts(jj,3,ii).^2];
%         colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii).^2,hopfB.pts(jj,3,ii).^2];
%         colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii).^2,hopfC.pts(jj,3,ii).^2];        
        colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii)*0.5+0.5,hopfA.pts(jj,3,ii)*0.5+0.5];
        colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii)*0.5+0.5,hopfB.pts(jj,3,ii)*0.5+0.5];
        colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii)*0.5+0.5,hopfC.pts(jj,3,ii)*0.5+0.5];
        plot3(hopfA.pts(jj,3,ii),hopfA.pts(jj,2,ii),hopfA.pts(jj,1,ii),'*','color',colorVA,'LineWidth',5)     
        plot3(hopfB.pts(jj,3,ii),hopfB.pts(jj,2,ii),hopfB.pts(jj,1,ii),'*','color',colorVB,'LineWidth',5)     
        plot3(hopfC.pts(jj,3,ii),hopfC.pts(jj,2,ii),hopfC.pts(jj,1,ii),'*','color',colorVC,'LineWidth',5)     
    end
    subplot(3,2,[1,2,3,4])
    axis equal;     grid off;    hold on;   axis off
    axis([-5 5 -5 5 -5 5])
    view([-30 15])
    for jj=1:nCircles
        if isnan(hopfA.pts(jj,1,ii)), continue; end
        colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii)*0.5+0.5,hopfA.pts(jj,3,ii)*0.5+0.5];
        colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii)*0.5+0.5,hopfB.pts(jj,3,ii)*0.5+0.5];
        colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii)*0.5+0.5,hopfC.pts(jj,3,ii)*0.5+0.5];
        plot3(hopfA.fibers(:,3,jj,ii),hopfA.fibers(:,2,jj,ii),hopfA.fibers(:,1,jj,ii),'-','color',hopfA.pts(jj,:,ii).^2,'LineWidth',2)
        plot3(hopfB.fibers(:,3,jj,ii),hopfB.fibers(:,2,jj,ii),hopfB.fibers(:,1,jj,ii),'-','color',hopfB.pts(jj,:,ii).^2,'LineWidth',2)
        plot3(hopfC.fibers(:,3,jj,ii),hopfC.fibers(:,2,jj,ii),hopfC.fibers(:,1,jj,ii),'-','color',hopfC.pts(jj,:,ii).^2,'LineWidth',2)
    end
    
end
