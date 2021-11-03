function rotor = rodrigues2(aXs,alpha)

%% Lie Algebra Matrices
LieZ = [0 -1 0;1 0 0;0 0 0];
LieY = [0 0 1;0 0 0;-1 0 0];
LieX = [0 0 0;0 0 -1;0 1 0];

%% Axis vector
Kvector = aXs(1)*LieX + aXs(2)*LieY + aXs(3)*LieZ;

%% Rotor
rotor = eye(3) + Kvector*sin(alpha) + Kvector*Kvector*(1-cos(alpha));
