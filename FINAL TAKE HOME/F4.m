clc;
clear all;
close all;
%% QUESTION 4
%%
reference = imread('Head-CT.png'); % original image
changing = imread('Head-CT-Transformed.png'); % image to be transformed
%% registration
[optimizer, metric] = imregconfig('monomodal');
transform = imregtform(changing, reference, 'rigid', optimizer, metric);
registered = imwarp(changing, transform, 'OutputView', imref2d(size(reference)));
%% difference calculation
difference = imabsdiff(reference, registered);
%% saving images
imwrite(registered, 'Head-CT-Registered.png');
imwrite(difference, 'Head-CT-Difference.png');
%% extracting the transformation parameters
rotation = rad2deg(asin(transform.T(2,1))); % Convert radians to degrees
translation_x = transform.T(3,1);
translation_y = transform.T(3,2);
fprintf('Optimum parameters for accurate transformation:\nRotation: %.2f degrees\n', rotation);
fprintf('Translation X: %.2f pixels\n', translation_x);
fprintf('Translation Y: %.2f pixels\n', translation_y);
%% applying transformation
%% figure 9
% rotation = 10;
% translation_x = -80;
% translation_y = 100;
%% figure 10
% rotation = 15;
% translation_x = -100;
% translation_y = 110;
%%
transform = affine2d([cosd(rotation) -sind(rotation) 0; sind(rotation) cosd(rotation) 0; translation_x translation_y 1]);
changing_registered = imwarp(changing, transform, 'OutputView', imref2d(size(reference)), 'SmoothEdges', true);
difference = imabsdiff(reference, changing_registered);
%% figuring results
figure(1)

subplot(141); 
imshow(reference);
title('REFERENCE IMAGE');

subplot(142); 
imshow(changing);
title('IMAGE TO BE TRANSFORMED');

subplot(143); 
imshow(changing_registered); 
title('REGISTERED TRANSFORMED IMAGE');

subplot(144)
imshow(difference);
title('IMAGE DIFFERENCE');
