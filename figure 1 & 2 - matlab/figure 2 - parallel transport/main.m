addpath('src')  
%% FIGURE 1
% unitary
[U, t] = get_unitary;
% bloch vector
[fig, eV, eB, TanVec, T, P] = bloch_frame(U, t, 60, -170);
make_figure1(fig, eV, eB, TanVec,'Figure #1')

%% FIGURE 2
%% The hopf fibration of the closed path
T=T(1:2^3:end); % polar angle
P=P(1:2^3:end); % aximuthal angle
Q = zeros(length(T), 4);
Q(:, 1) = cos(T/2).*cos(P/2);
Q(:, 2) = -sin(T/2).*sin(P/2);
Q(:, 3) = sin(T/2).*cos(P/2);
Q(:, 4) = cos(T/2).*sin(P/2);

%% FIGURE 2
% rotate the quaternion in the 6 Cayley axes
[Qli, Qlj, Qlk, Qri, Qrj, Qrk] = quaternion2circle(Q, 300);
% Hopf Projections - points on the 2 sphere:
[Rli, Rlj, Rlk, Rri, Rrj, Rrk] = quaternion2hopf(Qli, Qlj, Qlk, Qri, Qrj, Qrk);
%% Stereographic projection
fib_rk = stereo_project(Qrk);       % k - hopf projection
[sph.x,sph.y,sph.z] = sphere();     % 2sphere
figure(2);clf
make_figure2(fib_rk, Rrk(:,:,1), sph, 220)