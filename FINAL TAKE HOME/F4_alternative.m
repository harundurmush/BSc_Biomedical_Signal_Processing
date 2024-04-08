clc;
clear all;
close all;

fixed = imread('Head-CT.png');
moving = imread('Head-CT-Transformed.png'); 

if size(fixed, 3) == 3
    fixed = rgb2gray(fixed);
end
if size(moving, 3) == 3
    moving = rgb2gray(moving);
end

% Update rotation and translation ranges
rotations = -180:10:180; % -180 to 180 degrees rotation
translationsX = -180:10:180; % x-axis translation
translationsY = -180:10:180; % y-axis translation

bestResult = struct('Difference', inf, 'Registered', [], 'Rotation', 0, 'TranslationX', 0, 'TranslationY', 0);

% Find
for r = rotations
    for tx = translationsX
        for ty = translationsY
% Create the transformation matrix
            tform = affine2d([cosd(r) -sind(r) 0; ...
            sind(r) cosd(r) 0; ...
            tx ty 1]);
                    % Transform the moving image
        registered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)));
        
        % Calculate the difference image
        difference = imabsdiff(fixed, registered);
        differenceScore = sum(difference(:));
        
        % Update the best result if the current one is better
        if differenceScore < bestResult.Difference
            bestResult = struct('Difference', differenceScore, ...
                                'Registered', registered, ...
                                'Rotation', r, ...
                                'TranslationX', tx, ...
                                'TranslationY', ty);
        end
    end
end
end

figure('Name', 'Best Registration Result', 'NumberTitle', 'off');
subplot(1, 2, 1);
imshow(bestResult.Registered);
title(sprintf('Rotation: %d\nTx: %d, Ty: %d', ...
bestResult.Rotation, bestResult.TranslationX, bestResult.TranslationY));
subplot(1, 2, 2);
imshow(imabsdiff(fixed, bestResult.Registered), []);
title(sprintf('Difference: %d', bestResult.Difference));