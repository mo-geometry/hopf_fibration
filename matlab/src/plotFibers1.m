function plotFibers1(nCircles,stereo)

%% Plotting
% The figure !
[ sph.x, sph.y, sph.z ] = sphere();
figure(1);clf
subplot(3,2,[5,6])
axis equal;     grid off;    hold on;   axis off
% surf( sph.x,sph.y,sph.z, 'edgecolor', 'none' ) ; % eyeGlobe
surf( sph.x,sph.y,sph.z) ; % eyeGlobe
colormap winter
shading interp
for ii=1:nCircles
    plot3(stereo.pts(ii,3),stereo.pts(ii,2),stereo.pts(ii,1),'*','color',stereo.colorRange(ii,:),'LineWidth',5)    
end
view([200 30])
subplot(3,2,[1,2,3,4])
axis equal;     grid off;    hold on;   axis off
% plot the 6 projections 
for ii=1:nCircles
    plot3(stereo.r_hopf(:,3,ii),stereo.r_hopf(:,2,ii),stereo.r_hopf(:,1,ii),'-','color',stereo.colorRange(ii,:),'LineWidth',2)
end
view([270 10])