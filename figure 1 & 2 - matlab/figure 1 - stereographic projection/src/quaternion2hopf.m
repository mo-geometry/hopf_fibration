function [Rli, Rlj, Rlk, Rri, Rrj, Rrk] = quaternion2hopf(Qli, Qlj, Qlk, Qri, Qrj, Qrk)

% initialize 
a = size(Qli);
Rli = zeros(a(1),3,a(3));
Rlj = zeros(a(1),3,a(3));
Rlk = zeros(a(1),3,a(3));
Rri = zeros(a(1),3,a(3));
Rrj = zeros(a(1),3,a(3));
Rrk = zeros(a(1),3,a(3));

% Rli fill
Rli(:, 1, :) = Qli(:, 1, :).^2+Qli(:, 2, :).^2-Qli(:, 3, :).^2-Qli(:, 4, :).^2;
Rli(:, 2, :) = 2*(Qli(:, 2, :).*Qli(:, 3, :)-Qli(:, 1, :).*Qli(:, 4, :));
Rli(:, 3, :) = 2*(Qli(:, 2, :).*Qli(:, 4, :)+Qli(:, 1, :).*Qli(:, 3, :));
% Rlj fill
Rlj(:, 1, :) = 2*(Qlj(:, 2, :).*Qlj(:, 3, :)+Qlj(:, 1, :).*Qlj(:, 4, :));
Rlj(:, 2, :) = Qlj(:, 1, :).^2-Qlj(:, 2, :).^2+Qlj(:, 3, :).^2-Qlj(:, 4, :).^2;
Rlj(:, 3, :) = 2*(Qlj(:, 3, :).*Qlj(:, 4, :)-Qlj(:, 1, :).*Qlj(:, 2, :));
% Rlk fill
Rlk(:, 1, :) = 2*(Qlk(:, 2, :).*Qlk(:, 4, :)-Qlk(:, 1, :).*Qlk(:, 3, :));
Rlk(:, 2, :) = 2*(Qlk(:, 3, :).*Qlk(:, 4, :)+Qlk(:, 1, :).*Qlk(:, 2, :));
Rlk(:, 3, :) = Qlk(:, 1, :).^2-Qlk(:, 2, :).^2-Qlk(:, 3, :).^2+Qlk(:, 4, :).^2;
% Rri fill
Rri(:, 1, :) = Qri(:, 1, :).^2+Qri(:, 2, :).^2-Qri(:, 3, :).^2-Qri(:, 4, :).^2;
Rri(:, 2, :) = 2*(Qri(:, 2, :).*Qri(:, 3, :)+Qri(:, 1, :).*Qri(:, 4, :));
Rri(:, 3, :) = 2*(Qri(:, 2, :).*Qri(:, 4, :)-Qri(:, 1, :).*Qri(:, 3, :));
% Rrj fill
Rrj(:, 1, :) = 2*(Qrj(:, 2, :).*Qrj(:, 3, :)-Qrj(:, 1, :).*Qrj(:, 4, :));
Rrj(:, 2, :) = Qrj(:, 1, :).^2-Qrj(:, 2, :).^2+Qrj(:, 3, :).^2-Qrj(:, 4, :).^2;
Rrj(:, 3, :) = 2*(Qrj(:, 3, :).*Qrj(:, 4, :)+Qrj(:, 1, :).*Qrj(:, 2, :));
% Rrk fill
Rrk(:, 1, :) = 2*(Qrk(:, 2, :).*Qrk(:, 4, :)+Qrk(:, 1, :).*Qrk(:, 3, :));
Rrk(:, 2, :) = 2*(Qrk(:, 3, :).*Qrk(:, 4, :)-Qrk(:, 1, :).*Qrk(:, 2, :));
Rrk(:, 3, :) = Qrk(:, 1, :).^2-Qrk(:, 2, :).^2-Qrk(:, 3, :).^2+Qrk(:, 4, :).^2;
