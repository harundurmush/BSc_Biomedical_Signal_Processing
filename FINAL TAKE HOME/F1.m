clc;
clear all;
close all;
%% QUESTION 1
%%
I = imread("Head-CT.png"); % reading image file
if size(I,3) == 3
    graycoded_I = rgb2gray(I); % converting rgb to gray code
else
    graycoded_I = I;
end
%% step a
%% substep i, ii, iii
[G_x, G_y] = gradient(double(graycoded_I)); % obtaining gradient_x and gradient_y
abs_Gx = abs(G_x);
abs_Gy = abs(G_y);
gradient_I = abs_Gx + abs_Gy;
%% substep iv
filter_smooth = fspecial("average", [5 5]); % smooth filter
I_smooth = imfilter(graycoded_I, filter_smooth, "replicate");
[Gx_smooth, Gy_smooth] = gradient(double(I_smooth));
abs_Gx_smooth = abs(Gx_smooth);
abs_Gy_smooth = abs(Gy_smooth);
gradient_I_smooth = abs_Gx_smooth + abs_Gy_smooth;
%% step b
%%
I_marr = im2double(I);
% smoothening the image with a filter
filter_marr = [0 0 1 0 0; 0 1 2 1 0; 1 2 -16 2 1; 0 1 2 1 0; 0 0 1 0 0];
filtered_I = conv2(I_marr,filter_marr);
% finding the zero crossings
[row,column] = size(filtered_I);
zerocross = zeros([row,column]);
for i=2:row-1
    for j=2:column-1
        if (filtered_I(i,j)>0)
             if (filtered_I(i,j+1)>=0 && filtered_I(i,j-1)<0) || (filtered_I(i,j+1)<0 && filtered_I(i,j-1)>=0)
                             
                zerocross(i,j)= filtered_I(i,j+1);
                        
            elseif (filtered_I(i+1,j)>=0 && filtered_I(i-1,j)<0) || (filtered_I(i+1,j)<0 && filtered_I(i-1,j)>=0)
                    zerocross(i,j)= filtered_I(i,j+1);
            elseif (filtered_I(i+1,j+1)>=0 && filtered_I(i-1,j-1)<0) || (filtered_I(i+1,j+1)<0 && filtered_I(i-1,j-1)>=0)
                  zerocross(i,j)= filtered_I(i,j+1);
            elseif (filtered_I(i-1,j+1)>=0 && filtered_I(i+1,j-1)<0) || (filtered_I(i-1,j+1)<0 && filtered_I(i+1,j-1)>=0)
                  zerocross(i,j)=filtered_I(i,j+1);
            end
                        
        end
            
    end
end
output = im2uint8(zerocross);
% thresholding
output_thresholded = output > 10;
%% step c
%%
canny_edge = edge(graycoded_I, "canny");
%% figuring the results
%%
figure (1)

subplot(141)
imshow(graycoded_I, []);
title("ORIGINAL IMAGE");

subplot(142)
imshow(abs_Gx, []);
title("|G_X|");

subplot(143)
imshow(abs_Gy, []);
title("|G_Y|");

subplot(144)
imshow(gradient_I, []);
title("|G_X| + |G_Y|");

figure(2)

subplot(141)
imshow(graycoded_I, []);
title("ORIGINAL IMAGE");

subplot(142)
imshow(abs_Gx_smooth, []);
title("|G_X| AFTER SMOOTHING");

subplot(143)
imshow(abs_Gy_smooth, []);
title("|G_Y| AFTER SMOOTHING");

subplot(144)
imshow(gradient_I_smooth, []);
title("|G_X| + |G_Y| AFTER SMOOTHING");

figure(3)

subplot(131)
imshow(graycoded_I, []);
title("ORIGINAL IMAGE");

subplot(132)
imshow(output_thresholded, []);
title("MARR-HILDRETH");

subplot(133)
imshow(canny_edge, []);
title("CANNY");








