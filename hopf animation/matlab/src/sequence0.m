function hopf = sequence0(nCircles)

hopfA = fibersPts(getAngles0a(60));
hopfB = fibersPts(getAngles0b(nCircles,80));
hopfC = fibersPts(getAngles0c(nCircles,80));

[~,~,framesA] = size(hopfA.pts);
[~,~,framesB] = size(hopfB.pts);
[~,~,framesC] = size(hopfC.pts);

hopf.pts = nan(nCircles,3,framesA+framesB+framesC);
hopf.pts(1,:,1:framesA) = hopfA.pts;
hopf.pts(:,:,framesA+1:framesA+framesB) = hopfB.pts;
hopf.pts(:,:,framesA+framesB+1:framesA+framesB+framesC) = hopfC.pts;

hopf.Q = nan(nCircles,4,framesA+framesB+framesC);
hopf.Q(1,:,1:framesA) = hopfA.Q;
hopf.Q(:,:,framesA+1:framesA+framesB) = hopfB.Q;
hopf.Q(:,:,framesA+framesB+1:framesA+framesB+framesC) = hopfC.Q;


hopf.fibers = stereographicProjection(hopf.Q);

% hopf.stereo = nan(nCircles,3,framesA+framesB+framesC);
% hopf.stereo(1,:,1:framesA) = hopfA.stereo;
% hopf.stereo(:,:,framesA+1:framesA+framesB) = hopfB.stereo;
% hopf.stereo(:,:,framesA+framesB+1:framesA+framesB+framesC) = hopfC.stereo;