clc;
clear all;
close all;

fixed = imread("Head-CT.png");
moving = imread("Head-CT-Transformed.png");

rots = -20:10:20; % rotating angles to be tried
TX_values = -10:5:10; % translation in x-axis
TY_values = -10:5:10; % translation in y-axis

best_results = [];
num_bestresults = 5; % keeping best 5 results

threshold = 0.1; % Adjust the threshold as needed

for i = 1:length(rots)
    for j = 1:length(TX_values)
        for k = 1:length(TY_values)
            transform = affine2d([cosd(rots(i)) -sind(rots(i)) 0; ...
                              sind(rots(i)) cosd(rots(i)) 0; ...
                              TX_values(j) TY_values(k) 1]);
            
            registered = imwarp(moving, transform, 'OutputView', imref2d(size(fixed)));
            
            % Calculate the absolute difference between fixed and registered images
            difference = imabsdiff(fixed, registered);
            
            % Calculate the similarity score as the sum of squared differences
            similarity_score = sum(difference(:).^2);
            
            % Store the results
            best_results = [best_results; struct('Similarity', similarity_score, 'Registered', registered, 'Rotation', rots(i), 'TranslationX', TX_values(j), 'TranslationY', TY_values(k))];
        end
    end
end

% Sort the best results based on similarity
[~, order] = sort([best_results.Similarity]);
best_results = best_results(order);

% Create a figure for displaying results
figure('Name', 'Best Registration Results', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

for n = 1:num_bestresults
    % Create subplots in a grid
    subplot(2, num_bestresults, n);
    imshow(best_results(n).Registered);
    title(sprintf('Rotation: %d\nTx: %d, Ty: %d', ...
    best_results(n).Rotation, best_results(n).TranslationX, best_results(n).TranslationY));
    
    subplot(2, num_bestresults, num_bestresults + n);
    % Display the similarity score
    imshow(difference, []);
    title(sprintf('Similarity Score: %.4f', best_results(n).Similarity));
    
    % Check if similarity score is below the threshold
    if best_results(n).Similarity < threshold
        disp('Found a good registration:');
        disp(['Rotation: ', num2str(best_results(n).Rotation), ...
              ' Tx: ', num2str(best_results(n).TranslationX), ...
              ' Ty: ', num2str(best_results(n).TranslationY)]);
        break; % Stop when a good registration is found
    end
end
