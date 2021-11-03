function fibers = stereographicProjection(Q)

nPts = 300;
[nC,~,nF] = size(Q);

% Initialize the angle
lambda = linspace(0,2*pi,nPts)';

circleQ = zeros(nPts,4,nC,nF);

for ii=1:nF
    for jj=1:nC
        circleQ(:,1,jj,ii) = Q(jj,1,ii).*cos(lambda) + Q(jj,2,ii).*sin(lambda);
        circleQ(:,2,jj,ii) = Q(jj,2,ii).*cos(lambda) - Q(jj,1,ii).*sin(lambda);
        circleQ(:,3,jj,ii) = Q(jj,3,ii).*cos(lambda) - Q(jj,4,ii).*sin(lambda);
        circleQ(:,4,jj,ii) = Q(jj,4,ii).*cos(lambda) + Q(jj,3,ii).*sin(lambda);
    end
end

v1 = ones(size(circleQ(:,1,:,:)));
denom = v1-circleQ(:,1,:,:);
fibers = zeros(nPts,3,nC,nF);
fibers(:,1,:,:) = circleQ(:,2,:,:)./denom;
fibers(:,2,:,:) = circleQ(:,3,:,:)./denom;
fibers(:,3,:,:) = circleQ(:,4,:,:)./denom;

