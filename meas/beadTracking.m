function beadTracking(block)
%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

end

%% Coordinate system transfromation
%%  It is very important to note that it is assumed that the
%%  transfromation matrix H transforms the coordinates (x0,y0) from the
%%  full image from the CMOS chip while here we work with a cropped image
%%  and the coordinates (x1,y1) realted to the cropped image. The situation
%%  is illustrated below. That is why the parameters x_offset and y_offset
%%  are used in the functions below.
%%
%%   y0->--------------------------------------------
%%  x0    full image        |                        |
%%  |                       | <-x_offset             |
%%  |                       |                        |
%%  |            y1->--------------------            |
%%  |           x1                       |           | 
%%  |  y_offset |   cropped image        |           | 
%%  |<--------->|                        |           |
%%  |           |                        |           |  
%%  |            ------------------------            |  
%%   ------------------------------------------------

%% Function: elarray2implane ======================================
%% Abstract: Transforms coordintes from electrode array coordinate system
%%   to the image coordinate system.
%% Parameters:
%%   beadPos - [x, y, z] position of the beads (x, y, z are column vectors)
%%   H - transfromation matrix from image cooridnates (x0,y0) to electrode 
%%       array coordinates
%%   x_offset - x0 = x1 + x_offset
%%   y_offset - y0 = y1 + y_offset
function beadPos = elarray2implane(beadPos, H, x_offset, y_offset)
    if nargin < 3
        x_offset = 0;
        y_offset = 0;
    end
    
    tmp = H\[beadPos(:, 1)'; beadPos(:, 2)'; ones(1, size(beadPos, 1))];
    beadPos(:, 1) = (tmp(1,:)./tmp(3,:))' - x_offset;
    beadPos(:, 2) = (tmp(2,:)./tmp(3,:))' - y_offset;
end

%% Function: implane2elarray ======================================
%% Abstract: Transforms coordintes from image coordinate system to the
%%   electrode array coordinate system.
%% Parameters:
%%   beadPos - [x, y, z] position of the beads (x, y, z are column vectors)
%%   H - transfromation matrix from image cooridnates (x0,y0) to electrode 
%%       array coordinates
%%   x_offset - x0 = x1 + x_offset
%%   y_offset - y0 = y1 + y_offset
function beadPos = implane2elarray(beadPos, H, x_offset, y_offset)
    if nargin < 3
        x_offset = 0;
        y_offset = 0;
    end
    
    tmp = H*[beadPos(:, 1)' + x_offset; beadPos(:, 2)' + y_offset; ones(1, size(beadPos, 1))];
    beadPos(:, 1) = (tmp(1,:)./tmp(3,:))';
    beadPos(:, 2) = (tmp(2,:)./tmp(3,:))';
end

%% Function: genPropagationMatrix ======================================
%% Abstract:
%%   Generates propagation matrix for specified distance of back-propagation
%%
function [ Hq ] = genPropagationMatrix( img, z, dx, lambda, n)
    [M, N] = size(img);

    lM = dx*M;
    lN = dx*N;
    fx=-1/(2*dx):1/lN:1/(2*dx)-1/lN;
    fy=-1/(2*dx):1/lM:1/(2*dx)-1/lM;
    [FX, FY] = meshgrid(fx, fy);

    % Rayleigh-Sommerfeld propagator
    Hq = exp(1i*2*pi*z*n/lambda*sqrt(1-(lambda*FX/n).^2 -(lambda*FY/n).^2));
    Hq = Hq.*(sqrt(FX.^2+FY.^2) < (n/lambda));

    Hq = single(fftshift(Hq));

% end genPropagatorMatrix
end

%% Function: rsBackProp ==============================================
%% Abstract:
%%   Back-propagation of img1 by propagation matrix Hq1 and Hq2, and
%%      image img2 by propagation matrix Hq3
%%
function [Erz1, Erz2, Erz3] = rsBackProp( img1, img2, Hq1, Hq2, Hq3 )
% Implemented according to:
% [1] O. Mudanyali, D. Tseng, C. Oh, S. O. Isikman, I. Sencan, W. Bishara, C. Oztoprak, S. Seo, B. Khademhosseini, and A. Ozcan, “Compact, Light-weight and Cost-effective Microscope based on Lensless Incoherent Holography for Telemedicine Applications,” Lab Chip, vol. 10, no. 11, pp. 1417–1428, Jun. 2010.
% [2] S.-H. Lee and D. G. Grier, “Holographic microscopy of holographically trapped three-dimensional structures,” Optics express, vol. 15, no. 4, pp. 1505–1512, 2007.
% [3] F. C. Cheong, B. J. Krishnatreya, and D. G. Grier, “Strategies for three-dimensional particle tracking with holographic video microscopy,” Opt. Express, vol. 18, no. 13, pp. 13563–13573, 2010.
% Note: In [1] They use spatial frequency and in [2,3] spatial angular
% frequency
%
% The derivation of the formulae for the propagator (without refractive
% index n) can be found in (p. 60):
% [3]J. W. Goodman, Introduction to Fourier optics. McGraw-Hill, 1996.

Bq1 = fft2(img1);
Bq2 = fft2(img2);

Erz1 = Hq1.*Bq1;
Erz2 = Hq2.*Bq1;
Erz3 = Hq3.*Bq2;

Erz1 = ifft2(Erz1);
Erz2 = ifft2(Erz2);
Erz3 = ifft2(Erz3);

Erz1 = abs(Erz1);
Erz2 = abs(Erz2);
Erz3 = abs(Erz3);

% end rsBackPropMud
end

%% Function: findBead_normxcorr =========================================
%% Abstract:
%%   Locate beads in given image
%%
function [ beadPos ] = findBead_normxcorr( img, dx, beadDia, thresh )
% Normalization of the image
maxVal = max(img(:));
minVal = min(img(:));
img = uint8(double(img - minVal) / double(maxVal - minVal) * 255);

% Generate pattern
beadHalfWidth = ceil( beadDia / dx / 2  ) + 1;
beadPatternGen = 127*ones(2 * beadHalfWidth + 1, 2 * beadHalfWidth + 1);
[x, y] = meshgrid(-beadHalfWidth:beadHalfWidth, -beadHalfWidth:beadHalfWidth);
distMatrix = sqrt(x.^2 + y.^2);
beadPatternGen(distMatrix < ( beadDia / dx / 2 )) = 255;
beadPatternGen = uint8(beadPatternGen);

%% Normalized cross-correlation
% Generated pattern as template
crossCorrGen = normxcorr2(beadPatternGen, 255 - img);

% Search for local maximums
[xpeak, ypeak] = nonmaxsup2d(crossCorrGen, thresh);
xpeak = xpeak - beadHalfWidth;
ypeak = ypeak - beadHalfWidth;

beadPos = [ypeak, xpeak];
% end findBead_normxcorr
end

%% Function: nonmaxsup2d ===============================================
%% Abstract:
%%   Locate local maxiums greater than a specified threshold
%%
function [xpeak, ypeak] = nonmaxsup2d(response, thresh)
    H = size(response, 1);
    W = size(response, 2);

    xpeak = [];
    ypeak = [];
    
    for i = 2:(H-1)
        for j = 2:(W-1)            
            if (response(i, j) > thresh)                
                surPixels = response(i-1:i+1, j-1:j+1); 
                surPixels(2,2) = 0;
                if (response(i, j) > max(surPixels(:)))
                    xpeak = [xpeak; j];
                    ypeak = [ypeak; i];
                end
            end            
        end
    end    
% end nonmaxsup2d
end

%% Function: centerOfMass ===============================================
%% Abstract:
%%   Locate center of mass of an image. Only the pixels within an circle
%%   located at center of the image and touching the borders of the image
%%   are taken into calculation.
%%
function [ x, y ] = centerOfMass( img, power )

    x = single(0);
    y = single(0);
    total_intensity = single(0);
    w(1) = size(img, 1);
    w(2) = size(img, 2);
    R = max(w);
    for i=1:size(img, 1)
        for j=1:size(img, 2)
            if sqrt((i-w(1)/2-.5)^2+(j-w(2)/2-.5)^2) < R/2
                total_intensity = total_intensity + img(i, j).^power;
                x = x + i*img(i, j).^power;
                y = y + j*img(i, j).^power;
            end
        end
    end

    x = x / total_intensity;
    y = y / total_intensity;

% end centerOfMass
end


%% Function: trackBead ===============================================
%% Abstract:
%%   Locate beads according to their last position
%%
function beadPos = trackBead( img, beadPos, halfWindowSize )

    for l=1:size(beadPos, 1)
        for j = 1:numel(halfWindowSize)
            beadPosNew = beadPos(l, :);
            initilization = 1;
            k = 0;
            while (norm(beadPos(l, :) - beadPosNew) > 0.5 || initilization) && k < 10
                k = k + 1;
                initilization = 0;
                beadPos(l, :) = beadPosNew;
                
                hwWinSize = halfWindowSize(j);
                xMin = max(int32(beadPos(l, 1) - hwWinSize),       0);
                xMax = min(int32(beadPos(l, 1) + hwWinSize),        size(img, 1));
                yMin = max(int32(beadPos(l, 2) - hwWinSize),       0);
                yMax = min(int32(beadPos(l, 2) + hwWinSize),       size(img, 2));

                if xMin <= 0 || xMax > size(img, 1) || yMin <= 0 || yMax > size(img, 2) || ...
                        isnan(xMin) || isnan(xMax) || isnan(yMin) || isnan(yMax)
                    beadPos(l, 1) = NaN;
                    beadPos(l, 2) = NaN;
                    continue;
                end
                imgCropped = 1 - img(xMin:xMax, yMin:yMax);

                if var(imgCropped(:)) > 0.00001;
                    [x_cent, y_cent] = centerOfMass(imgCropped, 4);
                    
                    beadPosNew(1) = beadPos(l, 1) - hwWinSize - 1 + x_cent;
                    beadPosNew(2) = beadPos(l, 2) - hwWinSize - 1 + y_cent;
                else
                    disp('Bead lost!')
                    beadPos(l, 1) = NaN;
                    beadPos(l, 2) = NaN;
                    break;
                end
            end
        end
        beadPos(l, :) = beadPosNew;     
    end
% end findBead
end

%% Function: findBeadAngInit =============================================
%% Abstract:
%%   Locate beads in the image from oblique illumination. Called in the
%%   initialisation step.
%%
function beadPos = findBeadAngInit( img, beadPos, distRange, alpha, H, x_offset, y_offset)
%FINDBEADSTRAIGHT Summary of this function goes here
%   Detailed explanation goes here

for l=1:size(beadPos, 1)
    % Calculate the end points of line on which the center of the bead can
    % located
    cutoutXmin = beadPos(l,1) + distRange(1)*cos(alpha);
    cutoutXmax = beadPos(l,1) + distRange(2)*cos(alpha);
    cutoutYmin = beadPos(l,2) + distRange(1)*sin(alpha);
    cutoutYmax = beadPos(l,2) + distRange(2)*sin(alpha);
        
    % Convert the position from the el. array coordinate system to the
    % image coordinate system (oblique illumination cutout)
    cutoutPos = elarray2implane([cutoutXmin cutoutYmin; cutoutXmax cutoutYmax], H, x_offset, y_offset);
    
    cutoutXmin = cutoutPos(1,1);
    cutoutXmax = cutoutPos(2,1);
    cutoutYmin = cutoutPos(1,2);
    cutoutYmax = cutoutPos(2,2);
        
    % Find the bead
    imCutoutXmin = round(min(cutoutXmin, cutoutXmax));
    imCutoutXmax = round(max(cutoutXmin, cutoutXmax));
    imCutoutYmin = round(min(cutoutYmin, cutoutYmax));
    imCutoutYmax = round(max(cutoutYmin, cutoutYmax));
    
    % Check whether the coordinates are in correct bounds
    if ( imCutoutXmin < 0 || imCutoutXmax > size(img,1) || isnan(imCutoutXmin) || isnan(imCutoutXmax) ... %x coords
            || imCutoutYmin < 0 || imCutoutYmax > size(img,2)) || isnan(imCutoutYmin) || isnan(imCutoutYmax) % y coords || ...                 
        beadPos(l, 1) = NaN;
        beadPos(l, 2) = NaN;
        continue;
    end

    % Find 1 - long window
    imgCropped = 1 - img(imCutoutXmin:imCutoutXmax - 1, imCutoutYmin:imCutoutYmax - 1);    
    [x_cent, y_cent] = centerOfMass(imgCropped, 4);
    beadPos(l, 1) = x_cent + imCutoutXmin - 1;
    beadPos(l, 2) = y_cent + imCutoutYmin - 1;
        
    % Find 2 - smaller windowos
    beadPos(l, :) = trackBead(img, beadPos(l, :), [10 5]);
    
    % Convert the position from the the image coordinate system to el. arry coordinate system
    beadPos(l, :) = implane2elarray(beadPos(l, 1:2), H, x_offset, y_offset);
end
% end findBeadAngInit
end

%% Function: trackBeadsStr ===============================================
%% Abstract:
%%   Find beads in the image from the straight illumination. The positions
%%   is returned in the electrode array coordinate system.
%%
function beadPos = trackBeadsStr(img, beadPos, H, halfWinSize, x_offset, y_offset)
    % Convert the position from the el. array coordinate system to the
    % image coordinate system (straight illumination cutout)
    beadPos = elarray2implane(beadPos, H, x_offset, y_offset);

    % Update step of the tracking
    beadPos = trackBead(img, beadPos, halfWinSize);

    % Convert the position from the the image coordinate system to el. array coordinate system
    beadPos = implane2elarray(beadPos, H, x_offset, y_offset);
end

%% Function: trackBeadsAng ===============================================
%% Abstract:
%%   Find beads in the image from the oblique illumination. 
%% Arguments:
%%      img             - Image where the beads are tracked
%%      beadPosStr      - Actual position of the diffraction patterns from straight illumination in el. array coordinate system [um]
%%      prevHeight      - Previous estimate of the levitaiton height [um]
%%      halWinSize      - Size of the window around estimated bead position where the actual bead position is searched
%% .
%% .
function beadPos = trackBeadsAng(img, beadPosStr, prevHeights, halfWinSize, heightPolyCoeff, H, x_offset, y_offset, alpha)
    % Estimate the position if the diffraction pattern from oblique illumination in the el. array coordinate system
    dxy = prevHeights / heightPolyCoeff;
    beadPosAng_Est = [beadPosStr(:,1) + dxy*cos(alpha), beadPosStr(:,2) + dxy*sin(alpha)];

    % Convert the position from the el. array coordinate system to the
    % image coordinate system (oblique illumination cutout)
    beadPosAng_Est = elarray2implane(beadPosAng_Est, H, x_offset, y_offset);

    % Locate the bead diffraction patterns in the image cutout
    beadPos = trackBead(img, beadPosAng_Est, halfWinSize);

    % Convert the position from the the image coordinate system to el. arry coordinate system
    beadPos = implane2elarray(beadPos, H, x_offset, y_offset);
end

%% Function: estimateHeight ==============================================
%% Abstract:
%%   Estimate height of beads from position difference of their "shadows"
%%
function height = estimateHeight( beadPosStr, beadPosAng, heightPolyCoeff)
    % estimate height and plot it into the image
    height = heightPolyCoeff*sqrt(sum((beadPosStr(:,1:2) - beadPosAng(:,1:2)).^2,2));
%end estimateHeight
end

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C-Mex counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 2;
block.NumOutputPorts = 4;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Straight illumination image
block.InputPort(1).Dimensions        = [960 512];
block.InputPort(1).DatatypeID  = 1;  % single
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(1).SamplingMode = 'Sample';

% Oblique illumination image
block.InputPort(2).Dimensions        = [960 512];
block.InputPort(2).DatatypeID  = 1;  % single
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;
block.InputPort(2).SamplingMode = 'Sample';

% Reconstructed image
block.OutputPort(1).Dimensions       = [256 256];
block.OutputPort(1).DatatypeID  = 1; % single
block.OutputPort(1).Complexity  = 'Real';
block.OutputPort(1).SamplingMode = 'Sample';

% Bead Positions - el array coordinates (together with height estimate)
maxTrackedBeads     = block.DialogPrm(2).Data; % maximum number of tracked beads
block.OutputPort(2).Dimensions  = [maxTrackedBeads 3];
block.OutputPort(2).DatatypeID  = 0; % double
block.OutputPort(2).Complexity  = 'Real';
block.OutputPort(2).SamplingMode = 'Sample';

% Bead Positions - microscope-like image coordinates
block.OutputPort(3).Dimensions  = [maxTrackedBeads 2];
block.OutputPort(3).DatatypeID  = 0; % double
block.OutputPort(3).Complexity  = 'Real';
block.OutputPort(3).SamplingMode = 'Sample';

% Electrodes Positions
block.OutputPort(4).Dimensions  = [8 10];
block.OutputPort(4).DatatypeID  = 0; % double
block.OutputPort(4).Complexity  = 'Real';
block.OutputPort(4).SamplingMode = 'Sample';

% % Bead Positions (oblique illumination) - microscope-like image coordinates
% % block.OutputPort(5).Dimensions  = [maxTrackedBeads 2];
% block.OutputPort(5).Dimensions  = 1;
% block.OutputPort(5).DatatypeID  = 0; % double
% block.OutputPort(5).Complexity  = 'Real';
% block.OutputPort(5).SamplingMode = 'Sample';

% Register parameters
block.NumDialogPrms     = 5;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [-1 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
% block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
% block.RegBlockMethod('Update', @Update);
% block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required
% block.RegBlockMethod('CheckParameters', @CheckParam);

%end setup
end

%%
%% CheckParam:
%%
% function CheckParam(block)
%     %    TODO
% end

%%
%% PostPropagationSetup:
%%   Functionality    : Setup work areas and state variables. Can
%%                      also register run-time methods here
%%   Required         : No
%%   C-Mex counterpart: mdlSetWorkWidths
%%
function DoPostPropSetup(block)
maxtrackedBeads     = block.DialogPrm(2).Data; % maximum number of tracked beads

block.NumDworks = 20;

block.Dwork(1).Name            = 'Hq0';  % Hq_microscopelike
block.Dwork(1).Dimensions      = 256*256;
block.Dwork(1).DatatypeID      = 1;      % single
block.Dwork(1).Complexity      = 'Complex'; % complex
% block.Dwork(1).UsedAsDiscState = true;

block.Dwork(2).Name            = 'Hq1';  % Hq_farStr
block.Dwork(2).Dimensions      = 256*256;
block.Dwork(2).DatatypeID      = 1;      % single
block.Dwork(2).Complexity      = 'Complex'; % complex

block.Dwork(3).Name            = 'Hq2';  % Hq_farAng
block.Dwork(3).Dimensions      = 512*512;
block.Dwork(3).DatatypeID      = 1;      % single
block.Dwork(3).Complexity      = 'Complex'; % complex

block.Dwork(4).Name            = 'firstRun';  % flag indicataing an iteration is the first one or not (whether the initializetion is to be done)
block.Dwork(4).Dimensions      = 1;
block.Dwork(4).DatatypeID      = 8;      % bool
block.Dwork(4).Complexity      = 'Real'; % real

block.Dwork(5).Name            = 'reconstr_img';  % Hq_farAng
block.Dwork(5).Dimensions      = 256*256;
block.Dwork(5).DatatypeID      = 1;      % single
block.Dwork(5).Complexity      = 'Real'; % real

block.Dwork(6).Name            = 'beadPos';  % Position of Beads - ([x1 y1 h1; x2 y2 h2; ...]')(:), in el array coordinate system, in microns
block.Dwork(6).Dimensions      = 3*maxtrackedBeads;
block.Dwork(6).DatatypeID      = 0;      % double
block.Dwork(6).Complexity      = 'Real'; % real


block.Dwork(7).Name            = 'H_str';  % Homography matrix (transformation matrix from image coordinate system to el. array coordinate system) - straight illumination
block.Dwork(7).Dimensions      = 3*3;
block.Dwork(7).DatatypeID      = 0;      % double
block.Dwork(7).Complexity      = 'Real'; % real

block.Dwork(8).Name            = 'H_ang';  % Homography matrix (transformation matrix from image coordinate system to el. array coordinate system) - straight illumination
block.Dwork(8).Dimensions      = 3*3;
block.Dwork(8).DatatypeID      = 0;      % double
block.Dwork(8).Complexity      = 'Real'; % real

block.Dwork(9).Name            = 'ElPos';  % Electrodes position
block.Dwork(9).Dimensions      = 8*5*2; % 8 electrodes, 5 corners, 2 coordinates per corner
block.Dwork(9).DatatypeID      = 0;      % double
block.Dwork(9).Complexity      = 'Real'; % real

%% Paramters
block.Dwork(10).Name            = 'maxTrackedBeads';  % maximum number of tracked beads
block.Dwork(10).Dimensions      = 1; 
block.Dwork(10).DatatypeID      = 0;      % double
block.Dwork(10).Complexity      = 'Real'; % real

block.Dwork(11).Name            = 'dx';  % pixel size of the image chip
block.Dwork(11).Dimensions      = 1; 
block.Dwork(11).DatatypeID      = 0;      % double
block.Dwork(11).Complexity      = 'Real'; % real

block.Dwork(12).Name            = 'distRange';  % range of possible distance between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
block.Dwork(12).Dimensions      = 2; 
block.Dwork(12).DatatypeID      = 0;      % double
block.Dwork(12).Complexity      = 'Real'; % real

block.Dwork(13).Name            = 'alpha';  % angle between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
block.Dwork(13).Dimensions      = 1; 
block.Dwork(13).DatatypeID      = 0;      % double
block.Dwork(13).Complexity      = 'Real'; % real

block.Dwork(14).Name            = 'heightPolyCoeff';  % coefficients of function determining height of and bead based on distance between its diffraciton patterns - h = heightPolyCoeff(2) + heightPolyCoeff(1)*dx
block.Dwork(14).Dimensions      = 1; 
block.Dwork(14).DatatypeID      = 0;      % double
block.Dwork(14).Complexity      = 'Real'; % real

block.Dwork(15).Name            = 'y_offset';  % horizontal shift the 960*512 window
block.Dwork(15).Dimensions      = 1; 
block.Dwork(15).DatatypeID      = 0;      % double
block.Dwork(15).Complexity      = 'Real'; % real

block.Dwork(16).Name            = 'x_offset';  % vertical shift of the window for miscroscope-like reconstruction, x_offset(1) - shift of image crop for straight illumination, x_offset(2) - shift of image crop for oblique illumination
block.Dwork(16).Dimensions      = 2; 
block.Dwork(16).DatatypeID      = 0;      % double
block.Dwork(16).Complexity      = 'Real'; % real

block.Dwork(17).Name            = 'beadSize';  % Average size of the beads [um]
block.Dwork(17).Dimensions      = 1; 
block.Dwork(17).DatatypeID      = 0;      % double
block.Dwork(17).Complexity      = 'Real'; % real

block.Dwork(18).Name            = 'normxcorrTh';  % Threshold for normalized cross-correlation
block.Dwork(18).Dimensions      = 1; 
block.Dwork(18).DatatypeID      = 0;      % double
block.Dwork(18).Complexity      = 'Real'; % real

block.Dwork(19).Name            = 'markBeads';  % 0 - automaticaly find beads to track, 1 - mark beads to track
block.Dwork(19).Dimensions      = 1; 
block.Dwork(19).DatatypeID      = 8;      % bool
block.Dwork(19).Complexity      = 'Real'; % real

block.Dwork(20).Name            = 'trackedBeads';  % actual number of tracked beads
block.Dwork(20).DatatypeID      = 0;      % double
block.Dwork(20).Dimensions      = 1; 
block.Dwork(20).Complexity      = 'Real'; % real

% block.Dwork(21).Name            = 'beadPosAng_px';  % Position of Beads - ([x1 y1 h1; x2 y2 h2; ...]')(:), in el array coordinate system, in microns
% % block.Dwork(21).Dimensions      = 2*maxtrackedBeads;
% block.Dwork(21).Dimensions      = 1;
% block.Dwork(21).DatatypeID      = 0;      % double
% block.Dwork(21).Complexity      = 'Real'; % real

end

%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C-MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)
% Load params - from file
calibFileName = block.DialogPrm(5).Data;
if ~exist(calibFileName, 'file')
  % File does not exist.
  error('Warning: file does not exist:\n%s', calibFileName);
end
load(calibFileName) % Load calibration paramaters

% Load params - from function argument
beadSize            =  block.DialogPrm(1).Data; % Average size of the beads [um]
maxTrackedBeads     = block.DialogPrm(2).Data; % maximum number of tracked beads
y_offset            = 2*round(block.DialogPrm(3).Data/2); % horizontal shift the 960*512 window (must be even)
markBeads           = block.DialogPrm(4).Data; % 0 - automaticaly find beads to track, 1 - mark beads to track

% generate the propagation matrices
Hq_microscopelike   = genPropagationMatrix(zeros(256, 256), z_microscopelike, 2*dx, lambdaStr, refIndex);
Hq_farStr           = genPropagationMatrix(zeros(256, 256), z_farStr, 2*dx, lambdaStr, refIndex);
Hq_farAng           = genPropagationMatrix(zeros(512, 512), z_farAng, dx, lambdaAng, refIndex);

% Calculate x_offsets
x_offset(1) = round((750 - H_str(1,2)*(y_offset/2 + 125) - H_str(1,3)) / H_str(1,1) - 128);
x_offset(1) = 2*round(x_offset(1)/2); % must be even
x_offset(2) = round((750 - H_ang(1,2)*(y_offset + 256) - H_ang(1,3)) / H_ang(1,1) - 256);
x_offset(2) = 2*round(x_offset(2)/2); % must be even

% Store the propagation matrices to Dwork vectors
block.Dwork(1).Data = Hq_microscopelike(:);
block.Dwork(2).Data = Hq_farStr(:);
block.Dwork(3).Data = Hq_farAng(:);

% Store homography matices
block.Dwork(7).Data = H_str(:);
block.Dwork(8).Data = H_ang(:);

% Store params in Dwork vector
block.Dwork(10).Data = maxTrackedBeads      ; % maximum number of tracked beads
block.Dwork(11).Data = dx                   ; % pixel size of the image chip
block.Dwork(12).Data = distRange            ; % range of possible distance between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
block.Dwork(13).Data = alph                 ; % angle between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
block.Dwork(14).Data = heightPolyCoeff      ; % coefficients of function determining height of and bead based on distance between its diffraciton patterns - h = heightPolyCoeff*dx
block.Dwork(15).Data = y_offset             ; % horizontal shift the 960*512 window
block.Dwork(16).Data = x_offset             ; % vertical shift of the window for miscroscope-like reconstruction
block.Dwork(17).Data = beadSize             ; % Average size of the beads [um]
block.Dwork(18).Data = normxcorrTh          ; % Threshold for normalized cross-correlation
block.Dwork(19).Data = boolean(markBeads)   ; % 0 - automaticaly find beads to track, 1 - mark beads to track
block.Dwork(20).Data = maxTrackedBeads      ; % actual number of beads to track

% firstRun - set to TRUE
block.Dwork(4).Data = true;

% Generate polylines for electrode visualization
[x1, y1] = meshgrid(0:100:1500, [0, (6*250)]);
tmp = H_str\[x1(:)'; y1(:)'; ones(1, numel(x1))];
x2 = tmp(1,:)./tmp(3,:);
y2 = tmp(2,:)./tmp(3,:);
x2 = reshape(x2, size(x1))' - x_offset(1);
y2 = reshape(y2, size(y1))' - y_offset;

elPos = [];
for i=1:2:16
    % Fit a line going through the boundary points
    %   upper boundary of the electrode - y = k00*x + k01
    %   lower boundary of the electrode - y = k10*x + k11
    x00 = x2(i,1); x01 = x2(i,2);
    y00 = y2(i,1); y01 = y2(i,2);
    k00 = (y00 - y01) / (x00 - x01);
    k01 = (y01*x00 - y00*x01) / (x00 - x01);

    x10 = x2(i+1,1); x11 = x2(i+1,2);
    y10 = y2(i+1,1); y11 = y2(i+1,2);
    k10 = (y10 - y11) / (x10 - x11);
    k11 = (y11*x10 - y10*x11) / (x10 - x11);

    yUpLeft = 1;
    xUpLeft = (yUpLeft - k01)/k00;
    yUpRight = 255;
    xUpRight = (yUpRight - k01)/k00;

    yDwnLeft = 1;
    xDwnLeft = (yDwnLeft - k11)/k10;
    yDwnRight = 255;
    xDwnRight = (yDwnRight - k11)/k10;

    elPos = [elPos;
            yUpLeft, xUpLeft, yUpRight, xUpRight, yDwnRight, xDwnRight, yDwnLeft, xDwnLeft, yUpLeft, xUpLeft
            % y2(i,1), x2(i,1), y2(i,2), x2(i,2),       y2(i+1,2), x2(i+1,2),      y2(i+1,1), x2(i+1,1), y2(i,1), x2(i,1)
        ];
end
block.Dwork(9).Data = elPos(:);

% Save dialog params
save './expData/dialogParams.mat' 'x_offset' 'y_offset' 'maxTrackedBeads' 'beadSize'
save './expData/elPosPx.mat' 'elPos'


% block.Dwork(21).Data = 10;

%end InitializeConditions
end

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C-MEX counterpart: mdlOutputs
%%
function Outputs(block)
% tic
firstRun            = block.Dwork(4).Data;
H_str               = reshape(block.Dwork(7).Data, 3, 3); % Homography matrix (transformation matrix from image coordinate system to el. array coordinate system) - straight illumination
H_ang               = reshape(block.Dwork(8).Data, 3, 3);  % Homography matrix (transformation matrix from image coordinate system to el. array coordinate system) -  oblique illumination
maxTrackedBeads     = block.Dwork(10).Data; % maximum number of tracked beads
dx                  = block.Dwork(11).Data; % pixel size of the image chip
distRange           = block.Dwork(12).Data; % range of possible distance between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
alph                = block.Dwork(13).Data; % angle between diffraction pattern of an bead in the image from the straight illuminaiton and in the image from the oblique illumination
heightPolyCoeff     = block.Dwork(14).Data; % coefficients of function determining height of and bead based on distance between its diffraciton patterns - h = heightPolyCoeff*dx
y_offset            = block.Dwork(15).Data; % horizontal shift the 960*512 window
x_offset            = block.Dwork(16).Data; % horizontal shift the 960*512 window
beadSize            = block.Dwork(17).Data; % Average size of the beads [um]
normxcorrTh         = block.Dwork(18).Data; % Threshold for normalized cross-correlation
markBeads           = block.Dwork(19).Data; % markBeads
numOfBeadsToTrack   = block.Dwork(20).Data; % number of beads to track

% Reconstruct the captured images
img_green        = block.InputPort(1).Data((x_offset(2)):(x_offset(2) + 511), :);

img_red_raw = block.InputPort(2).Data(1:2:end, 2:2:end);
img_red = img_red_raw((x_offset(1)):(x_offset(1) + 255), :);

Hq_microscopelike   = reshape(block.Dwork(1).Data, 256, 256);
Hq_farStr           = reshape(block.Dwork(2).Data, 256, 256);
Hq_farAng           = reshape(block.Dwork(3).Data, 512, 512);

[img_microscopelike, img_farStr, img_farAng]  = rsBackProp(img_red, img_green, Hq_microscopelike, Hq_farStr, Hq_farAng);

if firstRun
    % firstRun - set to FALSE
    block.Dwork(4).Data = false;
    
    % Find bead in the numericaly reconstructed microscope-like image
    beadPos = single( findBead_normxcorr( img_microscopelike, 2*dx, beadSize, normxcorrTh ) );
    tmp = NaN*zeros(maxTrackedBeads, 3);
    
    if numel(beadPos) == 0
        disp('No beads were found!');
        % Store number of tracked beads
        block.Dwork(20).Data = 0;
    else
        % Refine the bead position estimates in
        beadPos = trackBead(img_farStr, beadPos, 7); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if markBeads
            fig = figure;
            imshow(img_microscopelike);

            r = 1.5*beadSize/4/dx;
            x_Circle = r*cos(0:pi/10:2*pi);
            y_Circle = r*sin(0:pi/10:2*pi);

            hold on;
            for i=1:size(beadPos,1)
                plot(beadPos(i,2)+x_Circle,  beadPos(i,1)+y_Circle, 'w')
            end
            hold off;

            beadPos = [];
            hold on;
            for i=1:maxTrackedBeads
                [x, y, but] = ginput(1);
                if but == 3
                    break;
                end
                plot(x, y, 'w*');
                beadPos = [beadPos; y, x];
            end
            hold off;            
            close(fig);
        end

        % Store number of beads to track
        block.Dwork(20).Data = size(beadPos, 1);

        % Convert the position from the the image coordinate system to el. array coordinate system
        beadPos = implane2elarray(beadPos, H_str, x_offset(1), y_offset/2);
                
        % find positions of the beads in image with oblique illumination
        beadPos_Ang = findBeadAngInit(img_farAng, beadPos, [-50 180], alph, H_ang, x_offset(2), y_offset);
        
        % Estimate heights of the beads
        height = estimateHeight(beadPos, beadPos_Ang, heightPolyCoeff);
        
        % Add height to the bead positions
        beadPos(:,3) = height;
        
        if size(beadPos,1) >= maxTrackedBeads
            tmp(1:maxTrackedBeads, :) = beadPos(1:maxTrackedBeads, :);
        else
            tmp(1:size(beadPos,1),:) = beadPos;
        end
    end
    
    tmp = tmp';
    block.Dwork(6).Data = tmp(:);
    
else
    beadPosPrev = reshape(block.Dwork(6).Data, 3, maxTrackedBeads)';
    beadPosPrev = beadPosPrev(1:numOfBeadsToTrack,:);
    
    % find positions of the beads in image with straight illumination
    beadPos = trackBeadsStr(img_farStr, beadPosPrev, H_str, [8 5], x_offset(1), y_offset/2);
    
    % find positions of the beads in image with oblique illumination
    beadPos_Ang = trackBeadsAng(img_farAng, beadPos, beadPosPrev(:,3), [10 5], heightPolyCoeff, H_ang, x_offset(2), y_offset, alph);
    
    % Estimate heights of the beads
    beadPos(:,3) = estimateHeight(beadPos, beadPos_Ang, heightPolyCoeff);
    
    beadPos = beadPos';
    block.Dwork(6).Data(1:3*numOfBeadsToTrack) = beadPos(:);

    % Bead position in image from oblique illumination
%     beadPos_Ang = elarray2implane(beadPos_Ang, H_ang, 0, 0);
%     block.Dwork(21).Data(1:2*numOfBeadsToTrack) = beadPos_Ang(:)';
end

% Zbytecny reshape - o par radek vyse se beadpos vektorizuje

beadPos = reshape(block.Dwork(6).Data, 3, maxTrackedBeads)';

beadPosInIm = elarray2implane(beadPos, H_str, x_offset(1), y_offset/2);

block.OutputPort(1).Data = img_microscopelike;
block.OutputPort(2).Data = beadPos;
block.OutputPort(3).Data = beadPosInIm(:, [2 1]);
block.OutputPort(4).Data = reshape(block.Dwork(9).Data, 8, 10);
% block.OutputPort(5).Data = reshape(block.Dwork(21).Data, 2, maxTrackedBeads)';
% block.OutputPort(5).Data = -1*block.Dwork(21).Data;
% block.Dwork(21).Data = -1*block.Dwork(21).Data;
% toc
%end Outputs
end

%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C-MEX counterpart: mdlTerminate
%%
function Terminate(block)

end
%end Terminate

