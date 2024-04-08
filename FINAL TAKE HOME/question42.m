clc;clear all;close all;
fixed = imread('Head-CT.png'); % The reference image
moving = imread('Head-CT-Transformed.png'); % The image to be registered

% Convert to grayscale if they are RGB
if size(fixed, 3) == 3
    fixed = rgb2gray(fixed);
end
if size(moving, 3) == 3
    moving = rgb2gray(moving);
end

% Initialize parameters for rotation (theta) and translation (tx, ty)
theta = 13; % Angle in degrees
tx = -100; % Translation in x
ty = 100; % Translation in y

% Adjust these parameters based on trial and error to find the best fit
% For example:
% theta = 10; 
% tx = 5;
% ty = -3;

% Perform transformation
tform = affine2d([cosd(theta) -sind(theta) 0; sind(theta) cosd(theta) 0; tx ty 1]);
movingRegistered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)), 'SmoothEdges', true);

% Calculate difference image
difference = imabsdiff(fixed, movingRegistered);

% Display images
figure;
subplot(1,3,1), imshow(fixed), title('Fixed Image');
subplot(1,3,2), imshow(movingRegistered), title('Registered Moving Image');
subplot(1,3,3), imshow(difference), title('Difference Image');
