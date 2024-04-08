%% Comments are added to the report (pdf file)
%%
clc;
clear all;
close all;
%% uploading images
%%
load('ecg_2.mat'); % ecg_hfn is from ecg_2.mat file
%% Question 2. Designing digital IIR filters using Butterworth and Chebychev filter approximations
%% step a) Butterworth and Chebyshev (Type I) filter parameters
%%
Fs = 1000; % sample frequency
fc = 30; % cutoff frequency
filter_order = 5; % choosing filter order
fc_normalized = fc / (Fs/2); % normalized cutoff frequency
pb_ripple = 2; % ripple in passband
[nom_butter, denom_butter] = butter(filter_order, fc_normalized, 'low'); % designing butterworth filter
ecg_filtered_butter = filter(nom_butter, denom_butter, ecg_hfn);
[nom_cheby1, denom_cheby1] = cheby1(filter_order, pb_ripple, fc_normalized, 'low'); % designing chebyshev type I filter
ecg_filtered_cheby1 = filter(nom_cheby1, denom_cheby1, ecg_hfn);
%% step a) time domain signals
%%
N_hfn = length(ecg_hfn);
Ts = 1/Fs;
t = 0:Ts:(N_hfn-1)*Ts;
figure (1)
subplot(211)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_filtered_butter);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Butterworth Approximation');
legend("original ecg signal","butterworth filtered ecg signal");
subplot(212)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_filtered_cheby1);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg signal","cheby (type I) filtered ecg signal");
%% step b) frequency domain signals
%%
ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_butter_filtered_freq = abs(fftshift(fft(ecg_filtered_butter,N_hfn)))/N_hfn;
ecg_cheby1_filtered_freq = abs(fftshift(fft(ecg_filtered_cheby1,N_hfn)))/N_hfn;
f = linspace(-Fs/2,Fs/2,N_hfn);
figure (2)
subplot(211)
plot(f, ecg_hfn_freq);
hold on;
plot(f, ecg_butter_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Butterworth Approximation');
legend("original ecg signal","butterworth filtered ecg signal");
subplot(212)
plot(f, ecg_hfn_freq);
hold on;
plot(f, ecg_cheby1_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type I) Approximation');
legend("original ecg signal","cheby (type I) filtered ecg signal");