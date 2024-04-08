clc;
clear all;
close all;
%% uploading images
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');
%% Question 1
%% step a
Fs = 1000;
Ts = 1/Fs;
%% step b
Td_hfn = length(ecg_hfn)/Fs;
Td_lfn = length(ecg_lfn)/Fs;
Td_noisy = length(ecg_noisy)/Fs;

%% step c

start_time = 1;  % seconds
end_time = 4;    % seconds
indices_range = (start_time * Fs)  : (end_time * Fs) - 1;
new_hfn1 = ecg_hfn(indices_range);
new_lfn1 = ecg_lfn(indices_range);
new_noisy1 = ecg_noisy(indices_range);
time1 = indices_range/Fs;

figure (1)
subplot(311)
plot(time1,new_hfn1);
title('From 2nd second to 5th second ecg-hfn(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-hfn(t)');

subplot(312)
plot(time1,new_lfn1);
title('From 2nd second to 5th second ecg-lfn(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-lfn(t)');

subplot(313)
plot(time1,new_noisy1);
title('From 2nd second to 5th second ecg-noisy(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-noisy(t)');
%% step d
% Typically, an Electrocardiogram (ECG) signal falls within the frequency band of 0.05 to 100 Hz. 
% This frequency range covers the relevant physiological information present in the electrical 
% activity of the heart.Here's a breakdown of the frequency bands in an ECG signal:
% Baseline Wander (0.05-0.5 Hz): This low-frequency component is associated with slow changes 
% in the ECG signal, such as variations in the baseline caused by respiration and body movement.
% Normal ECG Frequencies (0.5-50 Hz): The main frequency range of interest for ECG analysis is 
% from 0.5 to 50 Hz. This encompasses the P-waves, QRS complexes, and T-waves, representing the 
% electrical activity of the atria and ventricles of the heart.
% Powerline Interference (50 or 60 Hz): In environments with electrical interference, you might 
% observe powerline noise at 50 or 60 Hz, depending on the region.
%% step e
N = length(time1);
N_hfn1 = length(new_hfn1);
N_lfn1 = length(new_lfn1);
N_noisy1 = length(new_noisy1);
fvec_hfn1 = linspace(-Fs/2, Fs/2, N);
fvec_lfn1 = linspace(-Fs/2, Fs/2, N);
fvec_noisy1 = linspace(-Fs/2, Fs/2, N);
hfn_f1 = abs(fftshift(fft(new_hfn1,N_hfn1)))/N_hfn1;
lfn_f1 = abs(fftshift(fft(new_lfn1,N_lfn1)))/N_lfn1;
noisy_f1 = abs(fftshift(fft(new_noisy1,N_noisy1)))/N_noisy1;

figure(2)
subplot(311);
plot(fvec_hfn1,hfn_f1);
title('Magnitude Spectrum of the ecg-hfn(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-hfn(f)');

subplot(312);
plot(fvec_lfn1,lfn_f1);
title('Magnitude Spectrum of the ecg-lfn(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-lfn(f)');

subplot(313);
plot(fvec_noisy1,noisy_f1);
title('Magnitude Spectrum of the ecg-noisy(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-noisy(f)');

phase_hfn1 = unwrap(angle(fftshift(fft(new_hfn1,N))));
phase_lfn1 = unwrap(angle(fftshift(fft(new_lfn1,N))));
phase_noisy1 = unwrap(angle(fftshift(fft(new_noisy1,N))));

figure(3)
subplot(311);
plot(fvec_hfn1,phase_hfn1);
title('Phase Spectrum of the ecg-hfn(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-hfn(f)');

subplot(312);
plot(fvec_lfn1,phase_lfn1);
title('Phase Spectrum of the ecg-lfn(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-lfn(f)');

subplot(313);
plot(fvec_noisy1,phase_noisy1);
title('Phase Spectrum of the ecg-noisy(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-noisy(f)');
%% step f
start_time = 4;  % seconds
end_time = 7;    % seconds
indices_range2 = (start_time * Fs)  : (end_time * Fs) - 1;
new_hfn2 = ecg_hfn(indices_range2);
new_lfn2 = ecg_lfn(indices_range2);
new_noisy2 = ecg_noisy(indices_range2);
time2 = indices_range2/Fs;

figure (4)
subplot(311)
plot(time2,new_hfn2);
title('From 5th second to 8th second ecg-hfn(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-hfn(t)');

subplot(312)
plot(time2,new_lfn2);
title('From 5th second to 8th second ecg-lfn(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-lfn(t)');

subplot(313)
plot(time2,new_noisy2);
title('From 5th second to 8th second ecg-noisy(t)');
ylabel('volt (V)');
xlabel('time (s)');
legend('ecg-noisy(t)');

N_hfn2 = length(new_hfn2);
N_lfn2 = length(new_lfn2);
N_noisy2 = length(new_noisy2);
fvec_hfn2 = linspace(-Fs/2, Fs/2, N_hfn2);
fvec_lfn2 = linspace(-Fs/2, Fs/2, N_lfn2);
fvec_noisy2 = linspace(-Fs/2, Fs/2, N_noisy2);
hfn_f2 = abs(fftshift(fft(new_hfn2,N_hfn2)))/N_hfn2;
lfn_f2 = abs(fftshift(fft(new_lfn2,N_lfn2)))/N_lfn2;
noisy_f2 = abs(fftshift(fft(new_noisy2,N_noisy2)))/N_noisy2;

figure(5)
subplot(311);
plot(fvec_hfn2,hfn_f2);
title('Magnitude Spectrum of the ecg-hfn(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-hfn(f)');

subplot(312);
plot(fvec_lfn2,lfn_f2);
title('Magnitude Spectrum of the ecg-lfn(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-lfn(f)');

subplot(313);
plot(fvec_noisy2,noisy_f2);
title('Magnitude Spectrum of the ecg-noisy(f)');
ylabel('Amplitude');
xlabel('Frequency (Hz)');
legend('ecg-noisy(f)');

phase_hfn2 = unwrap(angle(fftshift(fft(new_hfn2,N_hfn2))));
phase_lfn2 = unwrap(angle(fftshift(fft(new_lfn2,N_lfn2))));
phase_noisy2 = unwrap(angle(fftshift(fft(new_noisy2,N_noisy2))));

figure(6)
subplot(311);
plot(fvec_hfn2,phase_hfn2);
title('Phase Spectrum of the ecg-hfn(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-hfn(f)');

subplot(312);
plot(fvec_lfn2,phase_lfn2);
title('Phase Spectrum of the ecg-lfn(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-lfn(f)');

subplot(313);
plot(fvec_noisy2,phase_noisy2);
title('Phase Spectrum of the ecg-noisy(f)');
ylabel('Angle');
xlabel('Frequency (Hz)');
legend('ecg-noisy(f)');
%% step g
%% Question 2
%%
% ecg_hfn will be used

% Sampling frequency
Fs = 1000;
Ts = 1/Fs;

% Design FIR filter
cutoff_freq = 50; % Adjust as needed
filter_order = 100; % Adjust as needed

% Design filter using window functions (e.g., Hamming and Blackman)
hamming_window = hamming(filter_order + 1);
blackman_window = blackman(filter_order + 1);

% Design FIR filters
fir_filter_hamming = fir1(filter_order, cutoff_freq/(Fs/2), 'low', hamming_window);
fir_filter_blackman = fir1(filter_order, cutoff_freq/(Fs/2), 'low', blackman_window);

% Apply filters to the ECG signal
filtered_signal_hamming = filter(fir_filter_hamming, 1, ecg_hfn);
filtered_signal_blackman = filter(fir_filter_blackman, 1, ecg_hfn);

% Plot original and filtered signals
figure(7)
subplot(2,1,1);
plot((0:length(ecg_hfn)-1)/Fs, ecg_hfn, 'b', 'LineWidth', 1.5);
hold on;
plot((0:length(filtered_signal_hamming)-1)/Fs, filtered_signal_hamming, 'r', 'LineWidth', 1.5);
title('Original vs. Filtered Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original', 'Filtered (Hamming)');

subplot(2,1,2);
plot((0:length(ecg_hfn)-1)/Fs, ecg_hfn, 'b', 'LineWidth', 1.5);
hold on;
plot((0:length(filtered_signal_blackman)-1)/Fs, filtered_signal_blackman, 'g', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original', 'Filtered (Blackman)');

% Plot frequency spectra
figure(8)
subplot(2,1,1);
plot_freq_spectrum(ecg_hfn, Fs, 'Original Signal');
subplot(2,1,2);
plot_freq_spectrum(filtered_signal_hamming, Fs, 'Filtered Signal (Hamming)');

% Function to plot frequency spectrum
function plot_freq_spectrum(signal, Fs, title_text)
    N = length(signal);
    fvec = linspace(-Fs/2, Fs/2, N);
    spectrum = abs(fftshift(fft(signal, N)))/N;
    plot(fvec, spectrum);
    title(title_text);
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
end





