function hopf = sequence3(nC,nF)

% Original
hopfA = fibersPts(getAngles3a(nC,nF));
hopfB = fibersPts(getAngles3b(nC,nF));
hopfC = fibersPts(getAngles3c(nC,3*nF));
% rotate
hopfB = fibersPts(angles4mPts(rotatePts3(hopfB.pts)));
% combine pts
hopf.pts(:,:,1:nF) = hopfA.pts;
hopf.pts(:,:,nF+1:2*nF) = hopfB.pts;
hopf.pts(:,:,2*nF+1:5*nF) = hopfC.pts;
% combine Quaternion
hopf.Q(:,:,1:nF) = hopfA.Q;
hopf.Q(:,:,nF+1:2*nF) = hopfB.Q;
hopf.Q(:,:,2*nF+1:5*nF) = hopfC.Q;

% stereographic projection
hopf.fibers = stereographicProjection(hopf.Q);

