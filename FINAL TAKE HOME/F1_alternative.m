clc;clear all;close all;
% Load the image
img = imread('Head-CT.png');

figure
imshow(img, []);
title('|gx| (Sobel x-direction)')
% Convert to grayscale if necessary
if size(img, 3) == 3
    img = rgb2gray(img);
end

% Define the Sobel masks
sobel_x = [-1 0 1; -2 0 2; -1 0 1];
sobel_y = [-1 -2 -1; 0 0 0; 1 2 1];

% Convolve the image with the Sobel masks
gx = conv2(double(img), sobel_x, 'same');
gy = conv2(double(img), sobel_y, 'same');

% Calculate the magnitudes
abs_gx = abs(gx);
abs_gy = abs(gy);

% Combine the gradients
combined_gradient = abs_gx + abs_gy;
figure
% Plot the magnitudes and combined gradient

imshow(abs_gx, []);
title('|gx| (Sobel x-direction)');
figure

imshow(abs_gy, []);
title('|gy| (Sobel y-direction)');


figure
imshow(combined_gradient, []);
title('|gx| + |gy| (Combined Sobel)');