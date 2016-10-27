function [ coord_elarray2implane, coord_implane2elarray, H ] = findElectrodes( img, markersInARow )

% Show the image and the message box asking for selecting the markers
figure;
imshow(img);
uiwait(msgbox(sprintf('Mark %d consecutive markers (rhombs) in each row form the left side to the right side.', markersInARow)));
mx = [];
my = [];
hold on
for i = 1:2*markersInARow
    [y, x] = ginput(1);
    
    mx = [mx; x];    
    my = [my; y];
    plot(y, x, 'b*')
    text(y, x, num2str(i))
end
hold off

save calibPoints mx my

%%

xi = mx;
yi = my;

% Se the corresponding coordinates of the markers in the el array coordinate system (im microns):
%  - (0, 0) is located below the left upmost marker on the upper edge of the electrode
%  - The markers are located 100+100*sqrt(2)/2 (~170.71) microns under and below electrodes (measured from the outter edge)
%  - Distance between two consecutive markers is 250 microns (center-to-center)
%  - Distance from the upper edge of the upmost electrode to the lower edge of the downmost electrode is 1500 microns (8 electrodes)
%  - x-axis is oriented downward and y-axis is oriented to the right
%  -> Thus, the markers in the upper row has coordinates ((-170.71, 0), (-170.71, 250), ...), in the lower row ((1670.71, 0), (1670.71, 250), ...)

xc = [-(100*sqrt(2)/2+100)*ones(markersInARow,1);
    1500 + (100*sqrt(2)/2+100)*ones(markersInARow,1)]; 

yc = 0:250:(markersInARow-1)*250;
yc = [yc(:); yc(:)];

imgSourad = [xi, yi];
elSourad = [xc, yc];

%% Zisserman's homography identification function
H_ziss = vgg_H_from_x_lin(imgSourad', elSourad')

% Show the electrodes in the image
[x1, y1] = meshgrid(0:100:1500, [0, ((markersInARow-1)*250)]);
tmp = H_ziss\[x1(:)'; y1(:)'; ones(1, numel(x1))];
x2 = tmp(1,:)./tmp(3,:);
y2 = tmp(2,:)./tmp(3,:);
x2 = reshape(x2, size(x1))';
y2 = reshape(y2, size(y1))';

hold on
for i=1:2:16
    plot( ...
    [y2(i,1), y2(i,2), y2(i+1,2), y2(i+1,1), y2(i,1)], ...
    [x2(i,1), x2(i,2), x2(i+1,2), x2(i+1,1), x2(i,1)], ...
    'y*-');
end
hold off

% Transform positions of the selected markers to the real positions
tmp = H_ziss*[xi(:)'; yi(:)'; ones(1, numel(xi))];
x2 = tmp(1,:)./tmp(3,:);
y2 = tmp(2,:)./tmp(3,:);

figure;
hold on;
plot(yc, xc ,'b*')
plot(y2, x2 ,'y*')
hold off

diffPos = sqrt((xc(:) - x2(:)).^2 + (yc(:) - y2(:)).^2);
sqrt(var(diffPos))

H = H_ziss;

    function [x2, y2] = elarray2implane(x1, y1)
        tmp = H\[x1(:)'; y1(:)'; ones(1, numel(x1))];
        x2 = tmp(1,:)./tmp(3,:);
        y2 = tmp(2,:)./tmp(3,:);
    end
    function [x2, y2] = implane2elarray(x1, y1)
        tmp = H*[x1(:); y1(:); ones(1, numel(x1))];
        x2 = tmp(1,:)./tmp(3,:);
        y2 = tmp(2,:)./tmp(3,:);
    end

coord_elarray2implane = @elarray2implane;
coord_implane2elarray = @implane2elarray;

end
