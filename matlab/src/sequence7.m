function [hopfA,hopfB,hopfC] = sequence7(nC,nF)

% Original
hopfA = fibersPts(getAngles7a(nC,nF));
hopfB = fibersPts(getAngles7b(nC,nF));
hopfC = fibersPts(getAngles7c(nC,nF));
% rotate
hopfA = fibersPts(angles4mPts(rotatePts7ac(hopfA.pts)));
hopfB = fibersPts(angles4mPts(rotatePts7b(hopfB.pts)));
hopfC = fibersPts(angles4mPts(rotatePts7ac(hopfC.pts)));
% stereographic projection
hopfA.fibers = stereographicProjection(hopfA.Q);
hopfB.fibers = stereographicProjection(hopfB.Q);
hopfC.fibers = stereographicProjection(hopfC.Q);

