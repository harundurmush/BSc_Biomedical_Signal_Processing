%% Comments are added to the report (pdf file)
%%
clc;
clear all;
close all;
%% uploading images
%%
load('ecg_2.mat'); % ecg_hfn is from ecg_2.mat file
%% Question 1. Applying a 10-point moving average filter
%% step a) time domain signals
%%
size_window = 10; % Defined window size of the filter
nom = ones(1,size_window); % nominator coefficient of the filter
denom = size_window; % denominator coefficient of the filter
ecg_hfn_filtered = filter(nom, denom, ecg_hfn);
N_hfn = length(ecg_hfn);
Fs = 1000; % sample frequency
Ts = 1/Fs;
t = 0:Ts:(N_hfn-1)*Ts;

figure (1)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_hfn_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg signal");
xlim([0 N_hfn*Ts]);
%% step b) frequency domain signals
%%
ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_hfn_filtered_freq = abs(fftshift(fft(ecg_hfn_filtered,N_hfn)))/N_hfn;
f = linspace(-Fs/2,Fs/2,N_hfn);

figure (2)
plot(f, ecg_hfn_freq);
hold on;
plot(f, ecg_hfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg signal");