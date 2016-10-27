%%
addpath('.\calibFunc')
close all;
clear all;
%%
lambdaStr         = 625e-9; % wavelength of green laser (the straight illumination) - 532
lambdaAng         = 525e-9; % wavelength of red laser (the oblique illumination) - 650
refIndex          = 1.44; % refractive index - it is chosen to be a constant, eventhough the light goes through different materials. As a result, the propagation heights DO NOT correspond to real heights and are chosen to give the most focused image.
dx                = 3.75e-6; % pixel size of the image chip

% Propagation distances
z_microscopelike    = 2400e-6; % Height where the numerically reconstructed image is in-focus
z_farStr            = 4500e-6; % Height where beads in the image from the straight illuminaiton (red) concentrate to one point
z_farAng            = 6100e-6; % Height where beads in the image from the oblique illuminaiton (green) concentrate to one point

normxcorrTh = 0.74;

save('calibParams.mat', 'lambdaStr', 'lambdaAng', 'dx', 'z_microscopelike', 'z_farStr', 'z_farAng', 'normxcorrTh', 'refIndex')

voltage = [ 10:-1:5 4.5:-0.5:2.5 0];
load('calibData01.mat')
for i = 1:numel(voltage)
    imwrite(img(:,:,:,i), sprintf('bottomView%03dV.png', 10*voltage(i)))
end

%% Measure heights of the bead form the side view images
height = measureHeightFromSideView_manual('./', voltage);
save heights_manual height

%% Find Electrodes
img = imread('bottomView100V.png');
% Find electrodes - Oblique
[ coord_elarray2implane_ang, coord_implane2elarray_ang, H_ang ] = findElectrodes( img(:,:,2), 6 );
waitforbuttonpress;
%%
% Find electrodes - Straight
imgStr = img(1:2:end, 2:2:end, 1);
[coord_elarray2implane_str, coord_implane2elarray_str, H_str ] = findElectrodes( imgStr, 6 );

save('calibParams.mat', 'H_str', 'H_ang', '-append')
waitforbuttonpress;

%% Differences in position of diffraction patterns from straight and oblique light sources
strWin1 = 8;
strWin2 = 5;
angWin1 = 10;
angWin2 = 4;
pow = 4;  
[posDif, alph, distRange] = measDistance_manual('./', voltage, H_str, H_ang, lambdaStr, lambdaAng, z_farStr, z_farAng, refIndex, dx, strWin1, strWin2, angWin1, angWin2, pow);

save('calibParams.mat', 'distRange', 'alph', '-append')

%%
heightPolyCoeff = fitPoly_height(height, posDif', voltage);

save('calibParams.mat', 'heightPolyCoeff', '-append')

%% Low pass filter
sysc = tf(1, [0.3 1]);
sysd = c2d(sysc, 0.1);

heightLPF.n = sysd.num{1};
heightLPF.d = sysd.den{1};

save('calibParams.mat', 'heightLPF', '-append')