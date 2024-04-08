clc;
clear all;
close all;
%% QUESTION 2
%%
I = imread('Head-CT.png');
graycoded_I = I;
level = graythresh(graycoded_I); % finding optimum threshold by Otsu's method
blacknwhite = imbinarize(graycoded_I, level);

figure(1)
imshow(blacknwhite);
title("SEGMENTATION OF THE IMAGE BY OTSU'S METHOD");

levels = multithresh(graycoded_I, 5);
segmented_I = imquantize(graycoded_I, levels);

figure(2)

subplot(131)
imshow(segmented_I == 1)
title("REGION 1");

subplot(132)
imshow(segmented_I == 2)
title("REGION 2");

subplot(133)
imshow(segmented_I == 3)
title("REGION 3");

I_smooth = imgaussfilt(graycoded_I, 5); % Gaussian filter with sigma value 5

level_smooth = graythresh(I_smooth);
blacknwhite_smooth = imbinarize(I_smooth, level_smooth);

figure(3)
imshow(blacknwhite_smooth);
title("SMOOTHED AND SEGMENTED IMAGE WITH OTSU'S METHOD");

levels_smooth = multithresh(I_smooth, 2);
segmented_I_smooth = imquantize(I_smooth, levels_smooth);

figure(4)

subplot(131)
imshow(segmented_I_smooth == 1)
title("SMOOTHED REGION 1");

subplot(132)
imshow(segmented_I_smooth == 2)
title("SMOOTHED REGION 2");

subplot(133)
imshow(segmented_I_smooth == 3)
title("SMOOTHED REGION 3");
