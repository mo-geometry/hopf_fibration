function hopf = sequence2(nC,nF)

% Original
hopf = fibersPts(getAngles2(nC,nF));
% rotate
hopf = fibersPts(angles4mPts(rotatePts2(hopf.pts)));
% stereographic projection
hopf.fibers = stereographicProjection(hopf.Q);
