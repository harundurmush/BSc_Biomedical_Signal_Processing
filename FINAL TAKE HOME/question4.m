clc;
clear all;
close all;
%% QUESTION 4
% Load images
fixed = imread('Head-CT.png'); % original image
moving = imread('Head-CT-Transformed.png'); % transformed image

% Perform registration
[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
registered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)));

% Calculate difference
difference = imabsdiff(fixed, registered);

% Save images
imwrite(registered, 'Head-CT-Registered.png');
imwrite(difference, 'Head-CT-Difference.png');

% Extract transformation parameters
% (assuming 'rigid' transformation which includes rotation and translation)
rotation = rad2deg(asin(tform.T(2,1))); % Convert radians to degrees
translation_x = tform.T(3,1);
translation_y = tform.T(3,2);

% Display the parameters
fprintf('Rotation: %.2f degrees\n', rotation);
fprintf('Translation X: %.2f pixels\n', translation_x);
fprintf('Translation Y: %.2f pixels\n', translation_y);

tform = affine2d([cosd(rotation) -sind(rotation) 0; sind(rotation) cosd(rotation) 0; translation_x translation_y 1]);
movingRegistered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)), 'SmoothEdges', true);

% Calculate difference image
difference = imabsdiff(fixed, movingRegistered);

% Display images
figure;
subplot(1,3,1), imshow(fixed), title('Fixed Image');
subplot(1,3,2), imshow(movingRegistered), title('Registered Moving Image');
