function ctr=plot2sphere3a(hopfA,hopfB,hopfC,frameString,str1)
ctr=0;
% initialize
[nCircles,~,nFrames] = size(hopfA.Q);
[ sph.x, sph.y, sph.z ] = sphere();
% generate sphere figure
sftX = 2;
sftY = -3;
sftZ = -3.5;

sX = sftX*ones(size(sph.x));
sY = sftY*ones(size(sph.x));
sZ = sftZ*ones(size(sph.x));

pX = sftX*ones(size(hopfA.pts(1,3,1)));
pY = sftY*ones(size(hopfA.pts(1,3,1)));
pZ = sftZ*ones(size(hopfA.pts(1,3,1)));
% plot points
for ii=1:nFrames
    
    
    
    %% plot
    h=figure(1);      clf;            
	set(gca,'color','w')
    hold on;
    axis equal;     grid off;    axis off
%     set(gca,'XTick',[]);set(gca,'YTick',[]) ;set(gca,'ZTick',[]) 
    surf( sph.x+sX,sph.y+sY,sph.z+sZ); 
    alpha(0.5)
    colormap gray
    shading interp
    view([30 15])
    for jj=1:nCircles
        if isnan(hopfA.pts(jj,1,ii)), continue; end
%         colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii).^2,hopfA.pts(jj,3,ii).^2];
%         colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii).^2,hopfB.pts(jj,3,ii).^2];
%         colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii).^2,hopfC.pts(jj,3,ii).^2];        
        colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii)*0.5+0.5,hopfA.pts(jj,3,ii)*0.5+0.5];
        colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii)*0.5+0.5,hopfB.pts(jj,3,ii)*0.5+0.5];
        colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii)*0.5+0.5,hopfC.pts(jj,3,ii)*0.5+0.5];
        plot3(hopfA.pts(jj,3,ii)+pX,hopfA.pts(jj,2,ii)+pY,hopfA.pts(jj,1,ii)+pZ,'*','color',colorVA,'LineWidth',3)     
        plot3(hopfB.pts(jj,3,ii)+pX,hopfB.pts(jj,2,ii)+pY,hopfB.pts(jj,1,ii)+pZ,'*','color',colorVB,'LineWidth',3)     
        plot3(hopfC.pts(jj,3,ii)+pX,hopfC.pts(jj,2,ii)+pY,hopfC.pts(jj,1,ii)+pZ,'*','color',colorVC,'LineWidth',3)     
    end
    for jj=1:nCircles
        if isnan(hopfA.pts(jj,1,ii)), continue; end
        colorVA = [hopfA.pts(jj,1,ii)*0.5+0.5,hopfA.pts(jj,2,ii)*0.5+0.5,hopfA.pts(jj,3,ii)*0.5+0.5];
        colorVB = [hopfB.pts(jj,1,ii)*0.5+0.5,hopfB.pts(jj,2,ii)*0.5+0.5,hopfB.pts(jj,3,ii)*0.5+0.5];
        colorVC = [hopfC.pts(jj,1,ii)*0.5+0.5,hopfC.pts(jj,2,ii)*0.5+0.5,hopfC.pts(jj,3,ii)*0.5+0.5];
        plot3(hopfA.fibers(:,3,jj,ii),hopfA.fibers(:,2,jj,ii),hopfA.fibers(:,1,jj,ii),'-','color',colorVA,'LineWidth',3)
        plot3(hopfB.fibers(:,3,jj,ii),hopfB.fibers(:,2,jj,ii),hopfB.fibers(:,1,jj,ii),'-','color',colorVB,'LineWidth',3)
        plot3(hopfC.fibers(:,3,jj,ii),hopfC.fibers(:,2,jj,ii),hopfC.fibers(:,1,jj,ii),'-','color',colorVC,'LineWidth',3)
    end
%     axis([-4.25 4.25 -4.25 4.25 -4.25 4.25])
    axis([-4.75 4.75 -4.75 4.75 -4.75 4.75])
    view([-30 15])
    
    ctr=ctr+1;
    if strcmp(str1,'writeVideo')
        saveas(h,[frameString,'/img',num2str(ctr,'%05.f'),'.png']);
    end
end
