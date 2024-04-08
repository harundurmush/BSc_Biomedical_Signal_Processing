clc;
clear all;
close all;
%% QUESTION 3
%%
I = imread("Head-CT.png");
if size(I,3) == 3
    graycoded_I = rgb2gray(I); % converting rgb to gray code
else
    graycoded_I = I;
end
reshaped_I = reshape(graycoded_I, [], 1);
k = [3, 9, 15]; % k values to be tried
figure;
for i = 1:length(k)
    [cluster_I_dx, ~] = kmeans(double(reshaped_I), k(i), "Distance", "sqEuclidean", "Replicates", 3);   
    clustered_I = reshape(cluster_I_dx, size(graycoded_I));
    max_value = max(clustered_I(:));
    clustered_I_scaled = uint8(clustered_I * (255 / max_value));
    
    subplot(1,3,i)
    imshow(clustered_I_scaled, []), title(["k-means clustering with k = ", num2str(k(i))]);
end
