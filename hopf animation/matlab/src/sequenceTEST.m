function hopf = sequenceTEST(nC)

% theta = pi/6*ones(nC,1);
phi0 = linspace(0,2*pi,nC)';
%% return angles
angles.theta(:,1,1) = pi/2*ones(nC,1);
angles.phi(:,1,1) = phi0;

%%
hopfA = fibersPts(angles);

[~,~,framesA] = size(hopfA.pts);

hopf.pts = nan(nC,3,framesA);
hopf.pts(:,:,1:framesA) = hopfA.pts;

hopf.Q = nan(nC,4,framesA);
hopf.Q(:,:,1:framesA) = hopfA.Q;

hopf.fibers = stereographicProjection(hopf.Q);