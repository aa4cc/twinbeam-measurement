maxTrackedBeads = 20;

% Lasers parameters
backProp.lambdaGreen = 532e-9; % wavelength of green laser (the straight illumination)
backProp.lambdaRed = 650e-9; % wavelength of red laser (the oblique illumination)
backProp.refIndex = 1.44; % refractive index - it is chosen to be a constant, eventhough the light goes through different materials. As a result, the propagation heights DO NOT correspond to real heights and are chosen to give the most focused image.
backProp.dx = 3.75e-6; % pixel size of the image chip

% Propagation distancex
backProp.z_microscopelike = 2400e-6; % Height where the numerically reconstructed image is in-focus
backProp.z_farStr = 4600e-6; % Height where beads in the image from the straight illuminaiton (green) concentrate to one point
backProp.z_farAng = 5400e-6; % Height where beads in the image from the oblique illuminaiton (red) concentrate to one point

% online searching params
threshStr = 0.15;
threshAng = 0.35;

% illumination params (distRange and alpha)
load('illuminationParams.mat');

% elevation height estimation params (fitted polynomial coefficients p - f(dist) = p(1) + p(2)*dist)
load('elheight_params.mat');
load('heightPlaneCoefs.mat')

% Camera & video recording initialization
y_offset = 250;
x_offset = 448;