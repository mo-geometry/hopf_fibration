function hopf = sequence1(nCircles,nFrames)

% Original
hopf1 = fibersPts(getAngles1(nCircles,nFrames));
% version a
hopfA = fibersPts(angles4mPts(rotatePts1a(hopf1.pts)));
hopfB = fibersPts(angles4mPts(rotatePts1b(hopf1.pts)));
hopfC = fibersPts(angles4mPts(rotatePts1c(hopf1.pts)));
% combine pts
hopf.pts(:,:,1:nFrames) = hopf1.pts;
hopf.pts(:,:,nFrames+1:2*nFrames) = hopfA.pts;
hopf.pts(:,:,2*nFrames+1:3*nFrames) = hopfB.pts;
hopf.pts(:,:,3*nFrames+1:4*nFrames) = hopfC.pts;
% combine Quaternion
hopf.Q(:,:,1:nFrames) = hopf1.Q;
hopf.Q(:,:,nFrames+1:2*nFrames) = hopfA.Q;
hopf.Q(:,:,2*nFrames+1:3*nFrames) = hopfB.Q;
hopf.Q(:,:,3*nFrames+1:4*nFrames) = hopfC.Q;
% stereographic projection
hopf.fibers = stereographicProjection(hopf.Q);









% hopfA = fibersPts(getAngles0a(120));
% hopfB = fibersPts(getAngles0b(nCircles,120));
% hopfC = fibersPts(getAngles0c(nCircles,80));
% 
% [~,~,framesA] = size(hopfA.pts);
% [~,~,framesB] = size(hopfB.pts);
% [~,~,framesC] = size(hopfC.pts);
% 
% hopf.pts = nan(nCircles,3,framesA+framesB+framesC);
% hopf.pts(1,:,1:framesA) = hopfA.pts;
% hopf.pts(:,:,framesA+1:framesA+framesB) = hopfB.pts;
% hopf.pts(:,:,framesA+framesB+1:framesA+framesB+framesC) = hopfC.pts;
% 
% hopf.Q = nan(nCircles,4,framesA+framesB+framesC);
% hopf.Q(1,:,1:framesA) = hopfA.Q;
% hopf.Q(:,:,framesA+1:framesA+framesB) = hopfB.Q;
% hopf.Q(:,:,framesA+framesB+1:framesA+framesB+framesC) = hopfC.Q;
% 
% 
% hopf.fibers = stereographicProjection(hopf.Q);
