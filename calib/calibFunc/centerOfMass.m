%% Function: centerOfMass ===============================================
%% Abstract:
%%   Locate center of mass in an image
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