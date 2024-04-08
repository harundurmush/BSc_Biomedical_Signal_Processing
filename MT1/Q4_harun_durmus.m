%% Comments are added to the report (pdf file)
%%
clc;
clear all;
close all;
%% uploading images
load('ecg_1.mat'); % ecg_lfn    is from ecg_1.mat file
load('ecg_2.mat'); % ecg_hfn    is from ecg_2.mat file
load('ecg_3.mat'); % ecg_noisy  is from ecg_3.mat file
%% Question 4. Power Spectral Density Analysis by using periodogram, modified periodogram, welch and levinson-durbin methods
%% step a) Using periodogram method (rectangular window)
fs = 1000; % sample frequency in Hz
N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);

[PSD_lfn_a, f_lfn_a] = periodogram(ecg_lfn, rectwin(N_lfn), N_lfn, fs);
[PSD_hfn_a, f_hfn_a] = periodogram(ecg_hfn, rectwin(N_hfn), (N_hfn+1)/2, fs);
[PSD_noisy_a, f_noisy_a] = periodogram(ecg_noisy, rectwin(N_noisy), N_noisy, fs);
PSD_lfn_dB_a = 10*log10(PSD_lfn_a);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_hfn_dB_a = 10*log10(PSD_hfn_a);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_noisy_dB_a = 10*log10(PSD_noisy_a); % converting magnitude of psd of the ecg_lfn signal to dB scale

figure(1)
subplot(311)
plot(f_lfn_a, PSD_lfn_dB_a);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Periodogram for estimating Power Spectral Density (PSD) of ecg-lfn');

subplot(312)
plot(f_hfn_a, PSD_hfn_dB_a);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Periodogram for estimating Power Spectral Density (PSD) of ecg-hfn');

subplot(313)
plot(f_noisy_a, PSD_noisy_dB_a);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Periodogram for estimating Power Spectral Density (PSD) of ecg-noisy');
%% step b) Using modified periodogram method (hamming window)
[PSD_lfn_b, f_lfn_b] = periodogram(ecg_lfn, hamming(N_lfn), N_lfn, fs);
[PSD_hfn_b, f_hfn_b] = periodogram(ecg_hfn, hamming(N_hfn), N_hfn, fs);
[PSD_noisy_b, f_noisy_b] = periodogram(ecg_noisy, hamming(N_noisy), N_noisy, fs);
PSD_lfn_dB_b = 10*log10(PSD_lfn_b);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_hfn_dB_b = 10*log10(PSD_hfn_b);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_noisy_dB_b = 10*log10(PSD_noisy_b); % converting magnitude of psd of the ecg_lfn signal to dB scale

figure(2)
subplot(311)
plot(f_lfn_b, PSD_lfn_dB_b);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Modified Periodogram for estimating Power Spectral Density (PSD) of ecg-lfn');

subplot(312)
plot(f_hfn_b, PSD_hfn_dB_b);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Modified Periodogram for estimating Power Spectral Density (PSD) of ecg-hfn');

subplot(313)
plot(f_noisy_b, PSD_noisy_dB_b);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title('Modified Periodogram for estimating Power Spectral Density (PSD) of ecg-noisy');
%% step c) Using welch's method
[PSD_lfn_c, f_lfn_c] = pwelch(ecg_lfn, hamming(256), 128, 256, fs);
[PSD_hfn_c, f_hfn_c] = pwelch(ecg_hfn, hamming(256), 128, 256, fs);
[PSD_noisy_c, f_noisy_c] = pwelch(ecg_noisy, hamming(256), 128, 256, fs);
PSD_lfn_dB_c = 10*log10(PSD_lfn_c);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_hfn_dB_c = 10*log10(PSD_hfn_c);     % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_noisy_dB_c = 10*log10(PSD_noisy_c); % converting magnitude of psd of the ecg_lfn signal to dB scale

figure(3)
subplot(311)
plot(f_lfn_c, PSD_lfn_dB_c);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title("Welch's method for estimating Power Spectral Density (PSD) of ecg-lfn");

subplot(312)
plot(f_hfn_c, PSD_hfn_dB_c);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title("Welch's method for estimating Power Spectral Density (PSD) of ecg-hfn");

subplot(313)
plot(f_noisy_c, PSD_noisy_dB_c);
grid on;
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');
title("Welch's method for estimating Power Spectral Density (PSD) of ecg-noisy");
%% step d) Using levin-durbin algorithm
order_lfn = 5;         % order of the AR model
order_hfn = 15;         % order of the AR model
order_noisy = 30;       % order of the AR model
nom_lfn = 1;        % nominator coefficient of freqz(.) is 1
nom_hfn = 1;        % nominator coefficient of freqz(.) is 1
nom_noisy = 1;      % nominator coefficient of freqz(.) is 1
[denom_lfn, error_lfn] = levinson(xcorr(ecg_lfn), order_lfn);
[denom_hfn, error_hfn] = levinson(xcorr(ecg_hfn), order_hfn);
[denom_noisy, error_noisy] = levinson(xcorr(ecg_noisy), order_noisy);
[PSD_lfn_d, f_lfn_d] = freqz(nom_lfn, denom_lfn, N_lfn, fs);
[PSD_hfn_d, f_hfn_d] = freqz(nom_hfn, denom_hfn, N_hfn, fs);
[PSD_noisy_d, f_noisy_d] = freqz(nom_noisy, denom_noisy, N_noisy, fs);
PSD_lfn_dB_d = 10*log10(abs(PSD_lfn_d).^2/error_lfn);          % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_hfn_dB_d = 10*log10(abs(PSD_hfn_d).^2/error_hfn);          % converting magnitude of psd of the ecg_lfn signal to dB scale
PSD_noisy_dB_d = 10*log10(abs(PSD_noisy_d).^2/error_noisy);    % converting magnitude of psd of the ecg_lfn signal to dB scale

figure (4)
subplot(311)
plot(f_lfn_d, PSD_lfn_dB_d);
grid on;
title('Levinson-Durbin Algorithm for estimating Power Spectral Density (PSD) of ecg-lfn');
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');

subplot(312)
plot(f_hfn_d, PSD_hfn_dB_d);
grid on;
title('Levinson-Durbin Algorithm for estimating Power Spectral Density (PSD) of ecg-hfn');
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');

subplot(313)
plot(f_noisy_d, PSD_noisy_dB_d);
grid on;
title('Levinson-Durbin Algorithm for estimating Power Spectral Density (PSD) of ecg-noisy');
xlabel('frequency (Hz)');
ylabel('power/frequency (dB/Hz)');