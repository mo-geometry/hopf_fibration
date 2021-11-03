function [U,t] = get_unitary
t = linspace(0,2*pi,2^12)';
% cayley matrices
cay.I = [1i 0;0 -1i];
cay.J = [0 1;-1 0];
cay.K = [0 1i;1i 0];
cay.e1 = [1 0;0 1];
% initialize
spinhalf = zeros(2,2,length(t));
mixed = zeros(2,2,length(t));
%% Define the unitary
for kk=1:length(t)
spinhalf(:,:,kk) = expm(cay.K*t(kk)/2)*expm(cay.J*t(kk)/2);
mixed(:,:,kk) = expm(-cay.I*t(kk))*expm(cay.K*t(kk)/2)*expm(cay.I*t(kk)/2)*expm(-cay.J*t(kk))*expm(-cay.I*t(kk));
end
% unitary
U=zeros(length(spinhalf),4);

for kk=1:length(spinhalf)
    U(kk,1) = real( spinhalf(1,1,kk) );
    U(kk,2) = imag( spinhalf(1,1,kk) );
    U(kk,3) = real( spinhalf(1,2,kk) );
    U(kk,4) = imag( spinhalf(1,2,kk) );    
end

