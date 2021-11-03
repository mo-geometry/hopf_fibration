function ctr=plot2sphere1a(hopf,folderString,str1)

% initialize
ctr=0;
[nCircles,~,nFrames] = size(hopf.Q);
[ sph.x, sph.y, sph.z ] = sphere();
% generate sphere figure
sftX = 2;
sftY = -3;
sftZ = -3.5;

sX = sftX*ones(size(sph.x));
sY = sftY*ones(size(sph.x));
sZ = sftZ*ones(size(sph.x));

pX = sftX*ones(size(hopf.pts(1,3,1)));
pY = sftY*ones(size(hopf.pts(1,3,1)));
pZ = sftZ*ones(size(hopf.pts(1,3,1)));
% plot points
for ii=1:nFrames
    
    
    
    %% plot
    h=figure(1);clf
    set(gca,'color','w')
    axis equal;     grid off;    hold on;   axis off
%     set(gca,'XTick',[]);set(gca,'YTick',[]) ;set(gca,'ZTick',[]) 
    surf( sph.x+sX,sph.y+sY,sph.z+sZ); 
    alpha(0.5)
    colormap gray
    shading interp
    view([30 15])
    for jj=1:nCircles
        if isnan(hopf.pts(jj,1,ii)), continue; end
%         colorV = [hopf.pts(jj,1,ii)*0.5+0.5,hopf.pts(jj,2,ii).^2,hopf.pts(jj,3,ii).^2];
        colorV = [hopf.pts(jj,1,ii)*0.5+0.5,hopf.pts(jj,2,ii)*0.5+0.5,hopf.pts(jj,3,ii)*0.5+0.5];
        plot3(hopf.pts(jj,3,ii)+pX,hopf.pts(jj,2,ii)+pY,hopf.pts(jj,1,ii)+pZ,'*','color',colorV,'LineWidth',3)     
    end
    for jj=1:nCircles
        if isnan(hopf.pts(jj,1,ii)), continue; end
        colorV = [hopf.pts(jj,1,ii)*0.5+0.5,hopf.pts(jj,2,ii)*0.5+0.5,hopf.pts(jj,3,ii)*0.5+0.5];
        plot3(hopf.fibers(:,3,jj,ii),hopf.fibers(:,2,jj,ii),hopf.fibers(:,1,jj,ii),'-','color',colorV,'LineWidth',3)
%         plot3(hopf.fibers(:,3,jj,ii),hopf.fibers(:,2,jj,ii),hopf.fibers(:,1,jj,ii),'r-','LineWidth',2)
    end
    axis([-4.25 4.25 -4.25 4.25 -4.25 4.25])
    view([-30 15])
    
        
    ctr=ctr+1;
    if strcmp(str1,'writeVideo')
        saveas(h,[ folderString,'/img',num2str(ctr,'%05.f'),'.png']);
    end
end
