function angles = angles4mPts(pts)

% pts = [cos(angles.theta),sin(angles.theta).*sin(angles.phi),sin(angles.theta).*cos(angles.phi)];
angles.theta = acos(pts(:,1,:));
angles.phi(:,1,:) = mod(atan2(pts(:,2,:),pts(:,3,:)),2*pi);