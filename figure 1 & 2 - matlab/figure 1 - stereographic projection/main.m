addpath('src') 

nPts = 20;
circle_res = 100;

for a=linspace(0.02,0.98, 50)
% declare omega
O = a*2*pi*ones(1, 3*nPts);
% declare phi
pp = linspace(0, 2*pi*0.85, nPts);
P = [pp, pp, pp];
% declare theta
T = pi*[0.33*ones(1, nPts), 0.5*ones(1, nPts), 0.66*ones(1, nPts)];
%% Quaternions
Q = zeros(3*nPts, 4);
Q(:, 1) = cos(T/2).*cos(P/2).*cos(O/2) - cos(T/2).*sin(P/2).*sin(O/2);
Q(:, 2) = sin(T/2).*cos(P/2).*sin(O/2) - sin(T/2).*sin(P/2).*cos(O/2);
Q(:, 3) = sin(T/2).*sin(P/2).*sin(O/2) + sin(T/2).*cos(P/2).*cos(O/2);
Q(:, 4) = cos(T/2).*sin(P/2).*cos(O/2) + cos(T/2).*cos(P/2).*sin(O/2);
%%
[Qli, Qlj, Qlk, Qri, Qrj, Qrk] = quaternion2circle(Q, circle_res);
% Hopf Projections (numerical check):
[Rli, Rlj, Rlk, Rri, Rrj, Rrk] = quaternion2hopf(Qli, Qlj, Qlk, Qri, Qrj, Qrk);
%% Stereographic projection
fib_li = stereo_project(Qli);
fib_lj = stereo_project(Qlj);
fib_lk = stereo_project(Qlk);
fib_ri = stereo_project(Qri);
fib_rj = stereo_project(Qrj);
fib_rk = stereo_project(Qrk);
%% Plot
[sph.x,sph.y,sph.z] = sphere();
h=figure(1);clf
subplot(3,2,1)
title('Left i hopf projection')
make_fig(fib_li, Rli(:,:,1), sph)
subplot(3,2,2)
title('Left j hopf projection')
make_fig(fib_lj, Rlj(:,:,1), sph)
subplot(3,2,3)
title('Left k hopf projection')
make_fig(fib_lk, Rlk(:,:,1), sph)
subplot(3,2,4)
title('Right i hopf projection')
make_fig(fib_ri, Rri(:,:,1), sph)
subplot(3,2,5)
title('Right j hopf projection')
make_fig(fib_rj, Rrj(:,:,1), sph)
subplot(3,2,6)
title('Right k hopf projection')
make_fig(fib_rk, Rrk(:,:,1), sph)
saveas(h,['stereo_hopf',num2str(round(a*100)),'.png'])
end