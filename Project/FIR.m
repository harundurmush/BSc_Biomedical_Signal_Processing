%% FIR filter (Moving Average Filter)
%%
clc;
clear all;
close all;
%% uploading image
load('ecg_1.mat'); % ecg_lfn is from ecg_1.mat file
load('ecg_2.mat'); % ecg_hfn is from ecg_2.mat file
load('ecg_3.mat'); % ecg_noisy is from ecg_3.mat file
%% Applying a 10-point moving average filter
%% time domain
size_window = 10; % Defined window size of the filter
nom = ones(1,size_window); % nominator coefficient of the filter
denom = size_window; % denominator coefficient of the filter
ecg_lfn_filtered = filter(nom, denom, ecg_lfn);
ecg_hfn_filtered = filter(nom, denom, ecg_hfn);
ecg_noisy_filtered = filter(nom, denom, ecg_noisy);
N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);
Fs = 1000; % sample frequency
Ts = 1/Fs;
t_lfn = 0:Ts:(N_lfn-1)*Ts;
t_hfn = 0:Ts:(N_hfn-1)*Ts;
t_noisy = 0:Ts:(N_noisy-1)*Ts;

figure (1)
subplot(311)
plot(t_lfn, ecg_lfn);
hold on;
plot(t_lfn, ecg_lfn_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg_lfn signal");
xlim([0 N_lfn*Ts]);

subplot(312)
plot(t_hfn, ecg_hfn);
hold on;
plot(t_hfn, ecg_hfn_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg_hfn signal");
xlim([0 N_hfn*Ts]);

subplot(313)
plot(t_noisy, ecg_noisy);
hold on;
plot(t_noisy, ecg_noisy_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg_noisy signal");
xlim([0 N_noisy*Ts]);
%% frequency domain
ecg_lfn_freq = abs(fftshift(fft(ecg_lfn,N_lfn)))/N_lfn;
ecg_lfn_filtered_freq = abs(fftshift(fft(ecg_lfn_filtered,N_lfn)))/N_lfn;
f_lfn = linspace(-Fs/2,Fs/2,N_lfn);

ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_hfn_filtered_freq = abs(fftshift(fft(ecg_hfn_filtered,N_hfn)))/N_hfn;
f_hfn = linspace(-Fs/2,Fs/2,N_hfn);

ecg_noisy_freq = abs(fftshift(fft(ecg_noisy,N_noisy)))/N_noisy;
ecg_noisy_filtered_freq = abs(fftshift(fft(ecg_noisy_filtered,N_noisy)))/N_noisy;
f_noisy = linspace(-Fs/2,Fs/2,N_noisy);

figure (2)
subplot(311)
plot(f_lfn, ecg_lfn_freq);
hold on;
plot(f_lfn, ecg_lfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg_lfn signal");

subplot(312)
plot(f_hfn, ecg_hfn_freq);
hold on;
plot(f_hfn, ecg_hfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg_hfn signal");

subplot(313)
plot(f_noisy, ecg_noisy_freq);
hold on;
plot(f_noisy, ecg_noisy_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg_noisy signal");