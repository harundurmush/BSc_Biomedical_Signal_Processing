%% Comments are added to the report (pdf file)
%%
clc;
clear all;
close all;
%% uploading images
load('ecg_1.mat'); % ecg_lfn    is from ecg_1.mat file
load('ecg_2.mat'); % ecg_hfn    is from ecg_2.mat file
load('ecg_3.mat'); % ecg_noisy  is from ecg_3.mat file
%% Question 3. About time-and-frequency analysis, Short-time Fourier Transform (STFT) and plotting the spectrogram.
%% step a) Using rectangular window
%%
n_window = 256;     % window length
n_overlap = 128;    % overlap length
Fs = 1000;          % sample frequency
[rect_lfn, f_lfn, t_lfn] = stft(ecg_lfn, Fs, 'Window', rectwin(n_window), 'OverlapLength', n_overlap);
[rect_hfn, f_hfn, t_hfn] = stft(ecg_hfn, Fs, 'Window', rectwin(n_window), 'OverlapLength', n_overlap);
[rect_noisy, f_noisy, t_noisy] = stft(ecg_noisy, Fs, 'Window', rectwin(n_window), 'OverlapLength', n_overlap);
% [s_lfn, f_lfn, t_lfn] = spectrogram(ecg_lfn, rectwin(n_window), n_overlap, [], Fs);
% [s_hfn, f_hfn, t_hfn] = spectrogram(ecg_hfn, rectwin(n_window), n_overlap, [], Fs);
% [s_noisy, f_noisy, t_noisy] = spectrogram(ecg_noisy, rectwin(n_window), n_overlap, [], Fs);
rect_lfn_dB = 10*log10(abs(rect_lfn));        % converting magnitude of STFT of the ecg_lfn signal to dB scale
rect_hfn_dB = 10*log10(abs(rect_hfn));        % converting magnitude of STFT of the ecg_hfn signal to dB scale
rect_noisy_dB = 10*log10(abs(rect_noisy));    % converting magnitude of STFT of the ecg_noisy signal to dB scale

figure(1)
subplot(311)
surf(t_lfn, f_lfn, rect_lfn_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-lfn with Rectangular Window');

subplot(312)
surf(t_hfn, f_hfn, rect_hfn_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-hfn with Rectangular Window');

subplot(313)
surf(t_noisy, f_noisy, rect_noisy_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-noisy with Rectangular Window');
%% step b) Using hamming window
%%
[hamming_lfn, f_lfn, t_lfn] = stft(ecg_lfn, Fs, 'Window', hamming(n_window), 'OverlapLength', n_overlap);
[hamming_hfn, f_hfn, t_hfn] = stft(ecg_hfn, Fs, 'Window', hamming(n_window), 'OverlapLength', n_overlap);
[hamming_noisy, f_noisy, t_noisy] = stft(ecg_noisy, Fs, 'Window', hamming(n_window), 'OverlapLength', n_overlap);
% [s_lfn, f_lfn, t_lfn] = spectrogram(ecg_lfn, rectwin(n_window), n_overlap, [], Fs);
% [s_hfn, f_hfn, t_hfn] = spectrogram(ecg_hfn, rectwin(n_window), n_overlap, [], Fs);
% [s_noisy, f_noisy, t_noisy] = spectrogram(ecg_noisy, rectwin(n_window), n_overlap, [], Fs);
hamming_lfn_dB = 10*log10(abs(hamming_lfn));        % converting magnitude of STFT of the ecg_lfn signal to dB scale
hamming_hfn_dB = 10*log10(abs(hamming_hfn));        % converting magnitude of STFT of the ecg_hfn signal to dB scale
hamming_noisy_dB = 10*log10(abs(hamming_noisy));    % converting magnitude of STFT of the ecg_noisy signal to dB scale

figure(2)
subplot(311)
surf(t_lfn, f_lfn, hamming_lfn_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-lfn with Hamming Window');

subplot(312)
surf(t_hfn, f_hfn, hamming_hfn_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-hfn with Hamming Window');

subplot(313)
surf(t_noisy, f_noisy, hamming_noisy_dB, 'EdgeColor', 'none');
axis xy; axis tight; colormap(jet); view(0, 90);
xlabel('time (s)');
ylabel('frequency (Hz)');
title('STFT of ecg-noisy with Hamming Window');