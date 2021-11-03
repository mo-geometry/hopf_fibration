function fibers = stereo_project(Q)

v1 = ones(size(Q(:,1,:)));
denom = v1-Q(:,1,:);
s_Q = size(Q);
fibers = zeros(s_Q(1),3,s_Q(3));
fibers(:,1,:) = Q(:,2,:)./denom;
fibers(:,2,:) = Q(:,3,:)./denom;
fibers(:,3,:) = Q(:,4,:)./denom;
