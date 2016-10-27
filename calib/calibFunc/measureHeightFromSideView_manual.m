function height = measureHeightFromSideView_manual(imgDir, voltages)
load('calibMatrix.mat');

fprintf('Mark the upper-left and lower-right corner of the Region of interest.\n\n')
imshow(imread(fullfile(imgDir, sprintf('sideView%03dV.png', voltages(1)*10))));
viewAreaBoundaries = ginput(2);
roi = [viewAreaBoundaries(1,1) viewAreaBoundaries(2,1) viewAreaBoundaries(1,2) viewAreaBoundaries(2,2)];

fprintf('Mark the bottom of the pool by marking one point on the bottom of the pool.\n\n')
imshow(imread(fullfile(imgDir, sprintf('sideView%03dV.png', voltages(end)*10))));
axis(roi)
[~, zeroHeightYcoord] = ginput(1);

height = zeros(size(voltages));
zeroHeightLinePoly = [0, zeroHeightYcoord];

fprintf('Mark the diameter of the bead. In other words, mark two points\nthat lie on a line going through the center of the bead and at\nthe same time lie on the boundary of the bead.\n\n')
[~, beadDia] = ginput(2);
beadDia = abs(beadDia(1)-beadDia(2));
xCirc = beadDia/2*cos(0:0.1:2*pi);
yCirc = beadDia/2*sin(0:0.1:2*pi);

%%
fprintf('Mark centers of the beads.\n')
fprintf('\t* left mouse button aproves the marked position\n\t* right mouse button disapproves the marker position and let\n\t  the user to mark the bead center again.\n')
%%

colormap gray;
for i=1:numel(voltages)
    img = imread(fullfile(imgDir, sprintf('sideView%03dV.png', voltages(i)*10)));
    imshow(img);
    axis(roi)
    
    while 1
        [x_pos, y_pos] = ginput(1);
        
        x_pos = mean(x_pos);
        y_pos = mean(y_pos);
        %     heightInPixels = size(img, 1) - y_pos;
        heightInPixels = polyval(zeroHeightLinePoly, x_pos) - y_pos;
        res = H*[x_pos x_pos; y_pos y_pos+heightInPixels; 1 1];
        height(i) = norm(res(1:2,1) - res(1:2,2)); %in microns
        
        hold on;
        plot([0 size(img,2)], polyval(zeroHeightLinePoly, [0 size(img,2)]), 'Color', [1 1 1]);
        plot(x_pos, y_pos, 'r*');
        plot([x_pos x_pos], [y_pos y_pos+heightInPixels], '-*y');
        text(x_pos, y_pos+heightInPixels/2, sprintf('%d um', round(height(i))), ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'right', ...
            'Color', [1 1 0]);
        expDesc(1) = {'Frequency 300 kHz'};
        expDesc(2) = {sprintf('Signal Amplitude %d V', i)};
        text(20, 20, expDesc, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'left', ...
            'Color', [1 1 0]);
        plot(x_pos + xCirc, y_pos + yCirc, 'y')
        hold off;
        
        [~, ~, but] = ginput(1);
        but
        if but ==1
            break;
        end
    end
    
    
        
    
end



