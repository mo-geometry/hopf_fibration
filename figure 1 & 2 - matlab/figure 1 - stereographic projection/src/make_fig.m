function make_fig(fibers, R, sph)
sft_x = -3; % move sphere x
sft_y = 4;  % move sphere y
sft_z = -1; % move sphere z
pts = length(R);
pt_colour = R/2+0.5;
%% subplot 1
axis equal;     grid on;    hold on;   axis on
% plot stereographic projection
for ii=1:pts
    plot3(squeeze(fibers(ii,1,:)),squeeze(fibers(ii,2,:)),squeeze(fibers(ii,3,:)),'-','color',pt_colour(ii,:),'LineWidth',4)
end
%% sphere
surf( sph.x+sft_x,sph.y+sft_y,sph.z+sft_z, 'FaceAlpha',0.5) ; % eyeGlobe
colormap winter
shading interp
% plot S2 points
for ii=1:pts
    plot3(R(ii,1)+sft_x,R(ii,2)+sft_y,R(ii,3)+sft_z,'*','color',pt_colour(ii,:),'LineWidth',3)    
end
% axis parameters
axis equal
xlim([-5,5])
ylim([-5,5])
zlim([-3,3])
xticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
xticklabels(['','','','','','','','','','',''])
yticks([-5,-4,-3,-2,-1,0,1,2,3,4,5])
yticklabels(['','','','','','','','','','',''])
zticks([-3,-2,-1,0,1,2,3])
zticklabels(['','','','','','',''])
xlabel('i')
ylabel('j')
zlabel('k')
view([140 10])