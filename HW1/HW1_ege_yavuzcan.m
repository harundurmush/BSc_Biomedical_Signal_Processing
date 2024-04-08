clc, clear all,close all;
%% QUESTION 1

load('ecg_1.mat'); 
load('ecg_2.mat');  
load('ecg_3.mat');  
ecg_1 = ecg_lfn;
ecg_2 = ecg_hfn;
ecg_3 = ecg_noisy;

fs = 1000;  % Hz

% a and b
samples_per_second = fs;


% b) Calculate the total recording times in seconds for each ECG data
total_time_ecg1 = length(ecg_1) / fs;
total_time_ecg2 = length(ecg_2) / fs;
total_time_ecg3 = length(ecg_3) / fs;

% c) Plot a total of 3 seconds of ECG data for each signal
start_time = 2;  % seconds
end_time = 5;    % seconds

% Calculate indices for the selected time range
indices_range = (start_time * fs)  : (end_time * fs) -1;

% Plot ECG data for the selected time range
figure;
subplot(3, 1, 1), plot(indices_range / fs, ecg_1(indices_range)), title('ECG 1');
subplot(3, 1, 2), plot(indices_range / fs, ecg_2(indices_range)), title('ECG 2');
subplot(3, 1, 3), plot(indices_range / fs, ecg_3(indices_range)), title('ECG 3');

% Calculate the number of samples for the 3-second fractions
samples_3_seconds = length(indices_range);

% d) Frequency band of typical ECG signal (around 0.05 to 150 Hz)

% e) Plot magnitude and phase spectrums for the selected time range
figure;
for i = 1:3
    % Calculate FFT
    fft_ecg = fft(ecg_1(indices_range));
    
    % Frequency axis
    freq_axis = linspace(0, fs, length(fft_ecg));

    % Plot magnitude spectrum
    subplot(3, 2, 2 * i - 1);
    plot(freq_axis, abs(fft_ecg));
    title(['Magnitude Spectrum - ECG ', num2str(i)]);
    
    % Plot phase spectrum
    subplot(3, 2, 2 * i);
    plot(freq_axis, angle(fft_ecg));
    title(['Phase Spectrum - ECG ', num2str(i)]);
end

% f) Repeat for the next 3 seconds
start_time = 5;  % seconds
end_time = 7;    % seconds

indices_range = (start_time * fs) + 1 : (end_time * fs);

% Calculate FFT for the new time range
fft_ecg_new_range = fft(ecg_1(indices_range));

% g) Compare magnitude and phase spectrums
figure;
subplot(2, 1, 1);
plot(freq_axis, abs(fft_ecg_new_range));
title('Magnitude Spectrum - New Time Range');

subplot(2, 1, 2);
plot(freq_axis, angle(fft_ecg_new_range));
title('Phase Spectrum - New Time Range');
