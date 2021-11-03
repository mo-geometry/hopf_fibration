function hopf = fibersPts(angles)

% points on the 2-sphere
hopf.pts = [cos(angles.theta),sin(angles.theta).*sin(angles.phi),sin(angles.theta).*cos(angles.phi)];
% quaternion
tmp1=sign(cos(angles.phi/2).*cos(angles.theta/2));
tmp2=[tmp1,tmp1,tmp1,tmp1];
hopf.Q = tmp2.*[cos(angles.phi/2).*cos(angles.theta/2),-sin(angles.phi/2).*cos(angles.theta/2),-cos(angles.phi/2).*sin(angles.theta/2),sin(angles.phi/2).*sin(angles.theta/2)];

% % stereographic projection
% v1 = ones(size(hopf.Q(:,1,:)));
% hopf.stereo = [hopf.Q(:,2,:)./(v1-hopf.Q(:,1,:)),hopf.Q(:,3,:)./(v1-hopf.Q(:,1,:)),hopf.Q(:,4,:)./(v1-hopf.Q(:,1,:))];