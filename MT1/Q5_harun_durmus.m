%% Comments are added to the report (pdf file)
%%
clc;
clear all;
close all;
%% uploading images
load('ecg_2.mat'); % ecg_hfn is from ecg_2.mat file
N_hfn = length(ecg_hfn);
fs = 1000; % sampling rate in Hz
ts = 1/fs;
t = 0:ts:(N_hfn-1)*ts;
%% Question 5. QRS Detection
%% low-pass Filter
% nominator part exist from the coefficients of x[n] and denominator part exist
% from the coefficients of y[n]:
% y[n] = 2y[n-1] - y[n-2] + x[n] - 2x[n-6] + x[n-12]
nom_low = [1,0,0,0,0,0,-2,0,0,0,0,0,1];
denom_low = [1,-2,1]; 
ecg_hfn_lowpass = filter(nom_low,denom_low,ecg_hfn);
%% high-pass filter
% y[n] = y[n-1] - x[n]/32 + x[n-16] - x[n-17] + x[n-32]/32
nom_high = [-1/32,zeros(1, 15),1,-1,zeros(1, 14),1/32]; 
denom_high = [1,-1];
ecg_hfn_highpass = filter(nom_high,denom_high,ecg_hfn_lowpass);
%% differentiator
% 8y[n] = 2x[n] + x[n-1] - x[n-3] - 2x[n-4]
nom_diff = [2,1,0,-1,-2]; 
denom_diff = 8;
ecg_hfn_diff = filter(nom_diff,denom_diff,ecg_hfn_highpass);
ecg_hfn_diff_0 = ecg_hfn_diff;
%% squaring
N_diff = length(ecg_hfn_diff);
for k=1:N_diff
    ecg_hfn_diff(k) = ecg_hfn_diff(k)^2;
end

figure (1)
subplot(221)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_hfn_lowpass);
xlabel("time (s)");
ylabel("amplitude");
subtitle("Signal Before and After Low-pass Filtering");
legend("before lowpass filter","after lowpass filter");
xlim([0 N_hfn*ts]);

subplot(222)
plot(t, ecg_hfn_lowpass);
hold on;
plot(t, ecg_hfn_highpass);
subtitle("Signal Before and After High-pass Filtering");
xlabel("time (s)");
ylabel("amplitude");
legend("before highpass filter","after highpass filter");
xlim([0 N_hfn*ts]);

subplot(223)
plot(t, ecg_hfn_highpass);
hold on;
plot(t, ecg_hfn_diff_0);
subtitle("Signal Before and After Differentiating");
xlabel("time (s)");
ylabel("amplitude");
legend("before differentiator","after differentiator");
xlim([0 N_hfn*ts]);

subplot(224)
plot(t, ecg_hfn_diff_0);
hold on;
plot(t, ecg_hfn_diff);
subtitle("Signal Before and After Squaring");
xlabel("time (s)");
ylabel("amplitude");
legend("before squaring","after squaring");
xlim([0 N_hfn*ts]);
%% moving average filter 
duration = 0.15; % 150 ms moving average filter parameter
size_window_1 = fs * duration; % defined window size of the filter
size_window_2 = 10; % the window size from the question 1
nom_maf_1 = ones(1,size_window_1); % nominator coefficient of the filter
nom_maf_2 = ones(1,size_window_2); % nominator coefficient of the filter
denom_maf_1 = size_window_1; % denominator coefficient of the filter
denom_maf_2 = size_window_2; % denominator coefficient of the filter
ecg_hfn_filtered_1 = filter(nom_maf_1, denom_maf_1, ecg_hfn_diff);
ecg_hfn_filtered_2 = filter(nom_maf_2, denom_maf_2, ecg_hfn_diff);
ecg_hfn_filtered_1 = ecg_hfn_filtered_1./2;
ecg_hfn_filtered_2 = ecg_hfn_filtered_2./20;

figure (2)
subplot(211)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_hfn_filtered_1);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title("Comparison of Two QRS Detected Signals With Different Window Sizes");
subtitle('Original and QRS Detected ECG Signal in time domain');
legend("original ecg signal","qrs detected ecg signal with window size 150");
xlim([0 N_hfn*ts]);

subplot(212)
plot(t, ecg_hfn);
hold on;
plot(t, ecg_hfn_filtered_2);
grid on;
xlabel("time (s)");
ylabel("amplitude");
subtitle('Original and QRS Detected ECG Signal in time domain');
legend("original ecg signal","qrs detected ecg signal with window size 10");
xlim([0 N_hfn*ts]);