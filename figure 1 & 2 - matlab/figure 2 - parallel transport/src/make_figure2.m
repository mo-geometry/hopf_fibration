function make_figure2(fibers, R, sph, a)

pts = length(R);
pt_colour = R/2+0.5;
%% subplot 1
% rotate the points on the 2 sphere 
a = a*pi/180;
rotor = [cos(a),sin(a),0; -sin(a), cos(a), 0; 0, 0, 1];
R = R*rotor;

%% subplot 2 
axis equal;     grid on;    hold on;   axis on
% plot the 6 projections 
for ii=1:pts
    plot3(squeeze(fibers(ii,1,:)),squeeze(fibers(ii,2,:)),squeeze(fibers(ii,3,:)),'-','color',pt_colour(ii,:),'LineWidth',4)
end

sft_x = -1.5;
sft_y = 2.0;
sft_z = -1.25;
a=5;

% 2sphere
surf( sph.x+sft_x,sph.y+sft_y,sph.z+sft_z, 'FaceAlpha',0.5) ; % eyeGlobe
colormap winter
shading interp
for ii=1:pts
    plot3(R(ii,1)+sft_x,R(ii,2)+sft_y,R(ii,3)+sft_z,'*','color',pt_colour(ii,:),'LineWidth',3)    
end

axis equal
xlim([-a,a])
ylim([-a,a])
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
