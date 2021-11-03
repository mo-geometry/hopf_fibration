function [hopfA,hopfB,hopfC] = sequence4(nC,nF)

% Original
hopfA = fibersPts(getAngles4a(nC,nF));
hopfB = fibersPts(getAngles4b(nC,nF));
hopfC = fibersPts(getAngles4c(nC,nF));
% % combine pts
% hopf.pts(:,:,1:nF) = hopfA.pts;
% hopf.pts(:,:,nF+1:2*nF) = hopfB.pts;
% hopf.pts(:,:,2*nF+1:3*nF) = hopfC.pts;
% % combine Quaternion
% hopf.Q(:,:,1:nF) = hopfA.Q;
% hopf.Q(:,:,nF+1:2*nF) = hopfB.Q;
% hopf.Q(:,:,2*nF+1:3*nF) = hopfC.Q;

% stereographic projection
hopfA.fibers = stereographicProjection(hopfA.Q);
hopfB.fibers = stereographicProjection(hopfB.Q);
hopfC.fibers = stereographicProjection(hopfC.Q);

