function [xM yM zM] = blochSp(n)

if nargin == 0,
    n = 4;
end

theta = linspace(0,pi,2^n)';
phi = linspace(0,2*pi,2^n)';

[th ph] = meshgrid(theta, phi);
xM = sin(th).*cos(ph);
yM = sin(th).*sin(ph);
zM = cos(th);