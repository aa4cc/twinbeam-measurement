function [posDif, ang, distRange] = measDistance_manual(imgDir, voltages, H_str, H_ang, lambdaStr, lambdaAng, z_farStr, z_farAng, n, dx, strWin1, strWin2, angWin1, angWin2, pow)

%% Back-propagation
% reconstruction distances
z_microscopelike = 2400e-6;

Hq_microscopelike = genPropagatorMatrix(zeros(960/2, 1280/2), z_microscopelike, 2*dx, lambdaStr, n);
Hq_farStr = genPropagatorMatrix(zeros(960/2, 1280/2), z_farStr, 2*dx, lambdaStr, n);
Hq_farAng = genPropagatorMatrix(zeros(960, 1280), z_farAng, dx, lambdaAng, n);

% load the image
img = imread(fullfile(imgDir, sprintf('bottomView%03dV.png', voltages(1)*10)));

% Crop and separate red and green components
imgGreen = img(:,:,2);
imgRed = img(1:2:end, 2:2:end, 1);

fprintf('Mark Regions of interests in both color channels containing the same bead.\n\n')
% Back-propagate
[~, img_farStr, img_farAng] = rsBackPropMud(imgRed, imgGreen, Hq_microscopelike, Hq_farStr, Hq_farAng);
imagesc(img_farStr);
str_axis = ginput(2);
imagesc(img_farAng);
ang_axis = ginput(2);

posDif = zeros(numel(voltages), 1);
ang = zeros(numel(voltages), 1);

fig = figure;
clf
% colormap gray
% colormap jet
colormap hsv
for i = 1:numel(voltages)
    % load the image
    img = imread(fullfile(imgDir, sprintf('bottomView%03dV.png', voltages(i)*10)));
    
    % Crop and separate red and green components
    imgGreen = img(:,:,2);
    imgRed = img(1:2:end, 2:2:end, 1);
    
    % Back-propagate
    [~, img_farStr, img_farAng] = rsBackPropMud(imgRed, imgGreen, Hq_microscopelike, Hq_farStr, Hq_farAng);
    
    % Manually find the center position of the bead
    imagesc(img_farStr)
    axis equal, axis tight
    axis([str_axis(1,1) str_axis(2,1) str_axis(1,2) str_axis(2,2)])
    title(sprintf('Straight, %2.1f V', voltages(i)));    
    hold on;
    [x, y] = ginput(1);
    beadPos = trackBead(img_farStr, [y x], [strWin1, strWin2], pow);
    x_str1(i) = beadPos(2);
    y_str1(i) = beadPos(1);
    plot(beadPos(2), beadPos(1), 'r*')
    hold off
    tmp = H_str*[y_str1(i); x_str1(i); 1];
    x_str = tmp(1,:)./tmp(3,:);
    y_str = tmp(2,:)./tmp(3,:);
    waitforbuttonpress
    
    imagesc(img_farAng)
    axis equal, axis tight
    axis([ang_axis(1,1) ang_axis(2,1) ang_axis(1,2) ang_axis(2,2)])
    title(sprintf('Oblique, %2.1f V', voltages(i)));
    hold on;
    [x, y] = ginput(1);
    beadPos = trackBead(img_farAng, [y x], [angWin1, angWin2], pow);
    x_ang1(i) = beadPos(2);
    y_ang1(i) = beadPos(1);
    plot(beadPos(2), beadPos(1), 'r*')
    hold off
    tmp = H_ang*[y_ang1(i); x_ang1(i); 1];
    x_ang = tmp(1,:)./tmp(3,:);
    y_ang = tmp(2,:)./tmp(3,:);
    waitforbuttonpress
    
    posDif(i) = norm([x_str-x_ang, y_str - y_ang]);
    ang(i) = mod(atan2(y_ang - y_str, x_ang - x_str), 2*pi);
    fprintf('Difference for voltage %2.1f: %f\n', voltages(i), posDif(i));
end

ang = mean(ang);
distRange = [0.5*min(posDif) 1.5*max(posDif)];

save beadPos x_str1 y_str1 x_ang1 y_ang1

close(fig)