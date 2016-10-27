function beadPos = trackBead( img, beadPos, halfWindowSize, pow )

for l=1:size(beadPos, 1)
    
%     figure(1)
%     imagesc(img)
%     axis equal
%     axis tight
%     axis([beadPos(2) - 4*halfWindowSize(1) beadPos(2) + 4*halfWindowSize(1) beadPos(1) - 4*halfWindowSize(2) beadPos(1) + 4*halfWindowSize(1)])
%     hold on
%     
%     rectangle('Position', [beadPos(l,2)-halfWindowSize(1), beadPos(l,1)-halfWindowSize(1), 2*halfWindowSize(1), 2*halfWindowSize(1)], 'EdgeColor', [0, 0, 0])
%     plot(beadPos(l,2), beadPos(l,1), '.k')
%     
%     c = [1 0 0; 0 1 0];
    for j = 1:numel(halfWindowSize)
        beadPosNew = beadPos(l, :);
        initilization = 1;
        k = 0;
        
        while (norm(beadPos(l, :) - beadPosNew) > 0.5 || initilization) && k <= 10
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
            imgCropped = 255 - img(xMin:xMax, yMin:yMax);
            
            if var(imgCropped(:)) > 0.00001;
                [x_cent, y_cent] = centerOfMass(imgCropped, pow);
                beadPosNew(1) = beadPos(l, 1) - hwWinSize - 1 + x_cent;
                beadPosNew(2) = beadPos(l, 2) - hwWinSize - 1 + y_cent;
                
%                 rectangle('Position', [beadPosNew(2)-halfWindowSize(j), beadPosNew(1)-halfWindowSize(j), 2*halfWindowSize(j), 2*halfWindowSize(j)], 'EdgeColor', (4-k+1)/4*c(j,:))
%                 plot(beadPosNew(2), beadPosNew(1), 'Marker', '.', 'Color', (4-k+1)/4*c(j,:))
            else
                disp('Bead lost!')
                beadPos(l, 1) = NaN;
                beadPos(l, 2) = NaN;
                break;
            end
        end
    end
    beadPos(l, :) = beadPosNew;
%     plot(beadPosNew(2), beadPosNew(1), 'yo')
%     hold off
end
end