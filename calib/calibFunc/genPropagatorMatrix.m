function [ Hq ] = genPropagatorMatrix( img, z, dx, lambda, n)
%GENPROPAGATORMATRIX Summary of this function goes here
%   Detailed explanation goes here

[M, N] = size(img);

lM = dx*M;
lN = dx*N;
fx=-1/(2*dx):1/lN:1/(2*dx)-1/lN;
fy=-1/(2*dx):1/lM:1/(2*dx)-1/lM;
[FX, FY] = meshgrid(fx, fy);

% Rayleigh-Sommerfeld
Hq = exp(1i*2*pi*z*n/lambda*sqrt(1-(lambda*FX/n).^2 -(lambda*FY/n).^2));
Hq = Hq.*(sqrt(FX.^2+FY.^2) < (n/lambda));

Hq = single(fftshift(Hq));

end

