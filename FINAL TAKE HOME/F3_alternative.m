clc;
clear all;
close all;

%%

img = imread('Head-CT.png');


grayImg = img;
if size(img, 3) == 3
    grayImg = rgb2gray(img);
end

reshapedImg = reshape(grayImg, [], 1);

% Farklı k değerleri ile k-means uygula
ks = [2, 6, 10]; % Denenecek k değerleri
for i = 1:length(ks)
    [clusterIdx, ~] = kmeans(double(reshapedImg), ks(i), 'Distance', 'sqEuclidean', 'Replicates', 3);
    
    clusteredImg = reshape(clusterIdx, size(grayImg));
    
  
    maxVal = max(clusteredImg(:));
    clusteredImgScaled = uint8(clusteredImg * (255 / maxVal));
    

    figure, 
    imshow(clusteredImgScaled, []), title(['k-means clustering with k = ', num2str(ks(i))]);
end


