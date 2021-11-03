function [Qli, Qlj, Qlk, Qri, Qrj, Qrk] = quaternion2circle(Q,n)

H = linspace(0,4*pi,n);
% initialize
Qli = zeros([size(Q),n]);
Qlj = zeros([size(Q),n]);
Qlk = zeros([size(Q),n]);
Qri = zeros([size(Q),n]);
Qrj = zeros([size(Q),n]);
Qrk = zeros([size(Q),n]);
% extract
a = Q(:, 1);
b = Q(:, 2);
c = Q(:, 3);
d = Q(:, 4);
% fill
shape = size(Qli);
for ii=1:shape(1)
% Qli fill
Qli(ii, 1, :) = a(ii) * cos(H/2) - b(ii) * sin(H/2);
Qli(ii, 2, :) = b(ii) * cos(H/2) + a(ii) * sin(H/2);
Qli(ii, 3, :) = c(ii) * cos(H/2) - d(ii) * sin(H/2);
Qli(ii, 4, :) = d(ii) * cos(H/2) + c(ii) * sin(H/2);
% Qlj fill
Qlj(ii, 1, :) = a(ii) * cos(H/2) - c(ii) * sin(H/2);
Qlj(ii, 2, :) = b(ii) * cos(H/2) + d(ii) * sin(H/2);
Qlj(ii, 3, :) = c(ii) * cos(H/2) + a(ii) * sin(H/2);
Qlj(ii, 4, :) = d(ii) * cos(H/2) - b(ii) * sin(H/2);
% Qlk fill
Qlk(ii, 1, :) = a(ii) * cos(H/2) - d(ii) * sin(H/2);
Qlk(ii, 2, :) = b(ii) * cos(H/2) - c(ii) * sin(H/2);
Qlk(ii, 3, :) = c(ii) * cos(H/2) + b(ii) * sin(H/2);
Qlk(ii, 4, :) = d(ii) * cos(H/2) + a(ii) * sin(H/2);
% Qri fill
Qri(ii, 1, :) = a(ii) * cos(H/2) + b(ii) * sin(H/2);
Qri(ii, 2, :) = b(ii) * cos(H/2) - a(ii) * sin(H/2);
Qri(ii, 3, :) = c(ii) * cos(H/2) - d(ii) * sin(H/2);
Qri(ii, 4, :) = d(ii) * cos(H/2) + c(ii) * sin(H/2);
% Qrj fill
Qrj(ii, 1, :) = a(ii) * cos(H/2) + c(ii) * sin(H/2);
Qrj(ii, 2, :) = b(ii) * cos(H/2) + d(ii) * sin(H/2);
Qrj(ii, 3, :) = c(ii) * cos(H/2) - a(ii) * sin(H/2);
Qrj(ii, 4, :) = d(ii) * cos(H/2) - b(ii) * sin(H/2);
% Qrk fill
Qrk(ii, 1, :) = a(ii) * cos(H/2) + d(ii) * sin(H/2);
Qrk(ii, 2, :) = b(ii) * cos(H/2) - c(ii) * sin(H/2);
Qrk(ii, 3, :) = c(ii) * cos(H/2) + b(ii) * sin(H/2);
Qrk(ii, 4, :) = d(ii) * cos(H/2) - a(ii) * sin(H/2);
end
