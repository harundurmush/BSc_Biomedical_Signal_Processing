%% Chebyshev Filter (Type I and Type II)
%%
clc;
clear all;
close all;
%% uploading image
load('ecg_1.mat'); % ecg_lfn is from ecg_1.mat file
load('ecg_2.mat'); % ecg_hfn is from ecg_2.mat file
load('ecg_3.mat'); % ecg_noisy is from ecg_3.mat file
%% Chebyshev Type I
%%
Fs = 1000; % sample frequency
fc = 30; % cutoff frequency
filter_order = 5; % choosing filter order
fc_normalized = fc / (Fs/2); % normalized cutoff frequency
pb_ripple = 2; % ripple in passband
[nom_cheby1, denom_cheby1] = cheby1(filter_order, pb_ripple, fc_normalized, 'low'); % designing chebyshev type I filter
ecg_filtered_cheby1 = filter(nom_cheby1, denom_cheby1, ecg_lfn);
ecg_filtered_cheby2 = filter(nom_cheby1, denom_cheby1, ecg_hfn);
ecg_filtered_cheby3 = filter(nom_cheby1, denom_cheby1, ecg_noisy);
%% time domain
N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);
Ts = 1/Fs;
t_lfn = 0:Ts:(N_lfn-1)*Ts;
t_hfn = 0:Ts:(N_hfn-1)*Ts;
t_noisy = 0:Ts:(N_noisy-1)*Ts;
figure;
subplot(311)
plot(t_lfn, ecg_lfn);
hold on;
plot(t_lfn, ecg_filtered_cheby1);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_lfn signal","cheby (type I) filtered ecg signal");

subplot(312)
plot(t_hfn, ecg_hfn);
hold on;
plot(t_hfn, ecg_filtered_cheby2);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_hfn signal","cheby (type I) filtered ecg signal");

subplot(313)
plot(t_noisy, ecg_noisy);
hold on;
plot(t_noisy, ecg_filtered_cheby3);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_noisy signal","cheby (type I) filtered ecg signal");

%% frequency domain
ecg_lfn_freq = abs(fftshift(fft(ecg_lfn,N_lfn)))/N_lfn;
ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_noisy_freq = abs(fftshift(fft(ecg_noisy,N_noisy)))/N_noisy;

ecg_cheby1_filtered_freq1 = abs(fftshift(fft(ecg_filtered_cheby1,N_lfn)))/N_lfn;
f_lfn = linspace(-Fs/2,Fs/2,N_lfn);
ecg_cheby1_filtered_freq2 = abs(fftshift(fft(ecg_filtered_cheby2,N_hfn)))/N_hfn;
f_hfn = linspace(-Fs/2,Fs/2,N_hfn);
ecg_cheby1_filtered_freq3 = abs(fftshift(fft(ecg_filtered_cheby3,N_noisy)))/N_noisy;
f_noisy = linspace(-Fs/2,Fs/2,N_noisy);

figure;
subplot(311)
plot(f_lfn, ecg_lfn_freq);
hold on;
plot(f_lfn, ecg_cheby1_filtered_freq1);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type I) Approximation');
legend("original ecg_lfn signal","cheby (type I) filtered ecg signal");

subplot(312)
plot(f_hfn, ecg_hfn_freq);
hold on;
plot(f_hfn, ecg_cheby1_filtered_freq2);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type I) Approximation');
legend("original ecg_hfn signal","cheby (type I) filtered ecg signal");

subplot(313)
plot(f_noisy, ecg_noisy_freq);
hold on;
plot(f_noisy, ecg_cheby1_filtered_freq3);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type I) Approximation');
legend("original ecg_noisy signal","cheby (type I) filtered ecg signal");
%% Chebyshev Type II
%%
Fs = 1000; % sample frequency
fc = 30; % cutoff frequency
filter_order = 5; % choosing filter order
fc_normalized = fc / (Fs/2); % normalized cutoff frequency
sb_ripple = 20; % ripple in stopband (dB), adjust as needed
% designing chebyshev type II filter
[nom_cheby2, denom_cheby2] = cheby2(filter_order, sb_ripple, fc_normalized, 'low');

% filtering the signals
ecg_filtered_cheby2_1 = filter(nom_cheby2, denom_cheby2, ecg_lfn);
ecg_filtered_cheby2_2 = filter(nom_cheby2, denom_cheby2, ecg_hfn);
ecg_filtered_cheby2_3 = filter(nom_cheby2, denom_cheby2, ecg_noisy);
%% time domain
N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);
Ts = 1/Fs;
t_lfn = 0:Ts:(N_lfn-1)*Ts;
t_hfn = 0:Ts:(N_hfn-1)*Ts;
t_noisy = 0:Ts:(N_noisy-1)*Ts;
figure;
subplot(311)
plot(t_lfn, ecg_lfn);
hold on;
plot(t_lfn, ecg_filtered_cheby2_1);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_lfn signal","cheby (type I) filtered ecg signal");

subplot(312)
plot(t_hfn, ecg_hfn);
hold on;
plot(t_hfn, ecg_filtered_cheby2_2);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_hfn signal","cheby (type I) filtered ecg signal");

subplot(313)
plot(t_noisy, ecg_noisy);
hold on;
plot(t_noisy, ecg_filtered_cheby2_3);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signal in time domain with Cheby (Type I) Approximation');
legend("original ecg_noisy signal","cheby (type I) filtered ecg signal");
%% frequency domain
ecg_lfn_freq = abs(fftshift(fft(ecg_lfn,N_lfn)))/N_lfn;
ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_noisy_freq = abs(fftshift(fft(ecg_noisy,N_noisy)))/N_noisy;

ecg_cheby2_filtered_freq1 = abs(fftshift(fft(ecg_filtered_cheby2_1,N_lfn)))/N_lfn;
f_lfn = linspace(-Fs/2,Fs/2,N_lfn);
ecg_cheby2_filtered_freq2 = abs(fftshift(fft(ecg_filtered_cheby2_2,N_hfn)))/N_hfn;
f_hfn = linspace(-Fs/2,Fs/2,N_hfn);
ecg_cheby2_filtered_freq3 = abs(fftshift(fft(ecg_filtered_cheby2_3,N_noisy)))/N_noisy;
f_noisy = linspace(-Fs/2,Fs/2,N_noisy);

figure;
subplot(311)
plot(f_lfn, ecg_lfn_freq);
hold on;
plot(f_lfn, ecg_cheby2_filtered_freq1);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type II) Approximation');
legend("original ecg_lfn signal","cheby (type II) filtered ecg signal");

subplot(312)
plot(f_hfn, ecg_hfn_freq);
hold on;
plot(f_hfn, ecg_cheby2_filtered_freq2);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type II) Approximation');
legend("original ecg_hfn signal","cheby (type II) filtered ecg signal");

subplot(313)
plot(f_noisy, ecg_noisy_freq);
hold on;
plot(f_noisy, ecg_cheby2_filtered_freq3);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain with Cheby (Type II) Approximation');
legend("original ecg_noisy signal","cheby (type II) filtered ecg signal");