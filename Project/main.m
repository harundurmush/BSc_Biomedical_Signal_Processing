%% **IT IS SUGGESTED THAT RUNNING THE CODE SECTION BY SECTION.**
%% LINE 9: FIR Filtering
%% LINE 115: IIR Filtering (Butterworth)
%% LINE 155: IIR Filtering (Chebyshev Type I)
%% LINE 257: IIR Filtering (Chebyshev Type  II)
%% LINE 363: Wavelet Transform
%% LINE 432: Wiener Filtering
%% LINE 486: Adaptive Filtering
%% FIR Filter
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

% Applying a 10-point moving average filter

% time domain

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

figure;
set(gcf,'color','w');
subplot(311)
plot(t_lfn, ecg_lfn);
hold on;
plot(t_lfn, ecg_lfn_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg-lfn signal");
xlim([0 N_lfn*Ts]);

subplot(312)
plot(t_hfn, ecg_hfn);
hold on;
plot(t_hfn, ecg_hfn_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg-hfn signal");
xlim([0 N_hfn*Ts]);

subplot(313)
plot(t_noisy, ecg_noisy);
hold on;
plot(t_noisy, ecg_noisy_filtered);
grid on;
xlabel("time (s)");
ylabel("amplitude");
title('Original and Filtered ECG Signals in time domain');
legend("original ecg signal","filtered ecg-noisy signal");
xlim([0 N_noisy*Ts]);

% frequency domain

ecg_lfn_freq = abs(fftshift(fft(ecg_lfn,N_lfn)))/N_lfn;
ecg_lfn_filtered_freq = abs(fftshift(fft(ecg_lfn_filtered,N_lfn)))/N_lfn;
f_lfn = linspace(-Fs/2,Fs/2,N_lfn);

ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N_hfn)))/N_hfn;
ecg_hfn_filtered_freq = abs(fftshift(fft(ecg_hfn_filtered,N_hfn)))/N_hfn;
f_hfn = linspace(-Fs/2,Fs/2,N_hfn);

ecg_noisy_freq = abs(fftshift(fft(ecg_noisy,N_noisy)))/N_noisy;
ecg_noisy_filtered_freq = abs(fftshift(fft(ecg_noisy_filtered,N_noisy)))/N_noisy;
f_noisy = linspace(-Fs/2,Fs/2,N_noisy);

figure;
set(gcf,'color','w');
subplot(311)
plot(f_lfn, ecg_lfn_freq);
hold on;
plot(f_lfn, ecg_lfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg-lfn signal");

subplot(312)
plot(f_hfn, ecg_hfn_freq);
hold on;
plot(f_hfn, ecg_hfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg-hfn signal");

subplot(313)
plot(f_noisy, ecg_noisy_freq);
hold on;
plot(f_noisy, ecg_noisy_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg-noisy signal");

%% IIR Filter (Butterworth)
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

Fs = 1000;
Ts = 1/Fs;
Fc = 50/Fs; %cutoff freq/Fs for normalization. cutoff freq determined by observing freq spectrum of the ecg signal
N = 4; %filter order
% by increasing order we get a sharper filter freq response supresses freqs above Fc better. 
% but increases phase shift and decreases amplitude more. it also requires more computational power so order of 4 is good enough
[b,a] = butter(N,Fc,"low");
ecg_butter_filtered = filter(b,a,ecg_hfn);

figure;
set(gcf,'color','w');
plot(ecg_hfn)
hold on
plot(ecg_butter_filtered)
title("ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_hfn","filtered ecg_hfn")

ecg_hfns=fftshift(abs(fft(ecg_hfn,length(ecg_hfn))/length(ecg_hfn)));
ecg_butter_filtereds=fftshift(abs(fft(ecg_butter_filtered,length(ecg_butter_filtered))/length(ecg_butter_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_butter_filtered));
f=transpose(f);

figure;
set(gcf,'color','w');
plot(f,ecg_hfns)
hold on
plot(f,ecg_butter_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_hfn","filtered ecg_hfn")

%% IIR Filter (Chebyshev Type I)
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

Fs = 1000; % sample frequency
fc = 30; % cutoff frequency
filter_order = 5; % choosing filter order
fc_normalized = fc / (Fs/2); % normalized cutoff frequency
pb_ripple = 2; % ripple in passband
[nom_cheby1, denom_cheby1] = cheby1(filter_order, pb_ripple, fc_normalized, 'low'); % designing chebyshev type I filter
ecg_filtered_cheby1 = filter(nom_cheby1, denom_cheby1, ecg_lfn);
ecg_filtered_cheby2 = filter(nom_cheby1, denom_cheby1, ecg_hfn);
ecg_filtered_cheby3 = filter(nom_cheby1, denom_cheby1, ecg_noisy);

% time domain

N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);
Ts = 1/Fs;
t_lfn = 0:Ts:(N_lfn-1)*Ts;
t_hfn = 0:Ts:(N_hfn-1)*Ts;
t_noisy = 0:Ts:(N_noisy-1)*Ts;

figure;
set(gcf,'color','w');
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

% frequency domain

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
set(gcf,'color','w');
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

%% IIR Filter (Chebyshev Type II)
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

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

% time domain

N_lfn = length(ecg_lfn);
N_hfn = length(ecg_hfn);
N_noisy = length(ecg_noisy);
Ts = 1/Fs;
t_lfn = 0:Ts:(N_lfn-1)*Ts;
t_hfn = 0:Ts:(N_hfn-1)*Ts;
t_noisy = 0:Ts:(N_noisy-1)*Ts;

figure;
set(gcf,'color','w');
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

% frequency domain

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
set(gcf,'color','w');
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

%% Denoising With Wavelet Transform (db4 or sym4)
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

Fs = 1000;
Ts = 1/Fs;
Fc = 50/Fs; %cutoff freq/Fs for normalization. cutoff freq determined by observing freq spectrum of the ecg signal
N = 4; %filter order
% by increasing order we get a sharper filter freq response supresses freqs above Fc better. 
% but increases phase shift and decreases amplitude more. it also requires more computational power so order of 4 is good enough
[b,a] = butter(N,Fc,"low");
ecg_butter_filtered = filter(b,a,ecg_hfn);

waveletType = 'sym4'; %'sym4'
noise_content = ecg_hfn - ecg_butter_filtered;
variance=var(noise_content);

%multilevel

level=5;
[c,l]=wavedec(ecg_hfn,level,waveletType);

j=1;
len=0;
for i=level:-1:1
    cD = detcoef(c,l, i); %detail coef in each level
    th=sqrt(2*variance*log(length(cD)));
    cD = wthresh(cD, 's', th);
    len=len+l(j);
    c(len+1:len+length(cD))=cD;
    j=j+1;
end

denoised_ecg_wavelet=waverec(c,l,waveletType);
denoised_ecg_wavelet_s=fftshift(abs(fft(denoised_ecg_wavelet,length(denoised_ecg_wavelet))/length(denoised_ecg_wavelet)));


% Plot original and denoised ECG signal

figure;
set(gcf,'color','w');
plot(ecg_hfn);
title('Original ECG Signal');
xlabel('Sample Number');
ylabel('Amplitude');
hold on
plot(denoised_ecg_wavelet);
title('Denoised ECG Signal using Wavelet (db4) with Custom Thresholding');
xlabel('Sample Number');
ylabel('Amplitude');
legend("ecg","denoised")

%freq domain plots
ecg_hfns=fftshift(abs(fft(ecg_hfn,length(ecg_hfn))/length(ecg_hfn)));
f=linspace(-Fs/2,Fs/2,length(ecg_butter_filtered));
f=transpose(f);

figure;
set(gcf,'color','w');
plot(f,ecg_hfns)
hold on
plot(f,denoised_ecg_wavelet_s)    
title("frequency spectrum of ecg signal and its filtered version with wavelet")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_hfn","filtered ecg_hfn")

%% Wiener Filter
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

Fs = 1000;
Ts = 1/Fs;
Fc = 50/Fs; %cutoff freq/Fs for normalization. cutoff freq determined by observing freq spectrum of the ecg signal
N = 4; %filter order
% by increasing order we get a sharper filter freq response supresses freqs above Fc better. 
% but increases phase shift and decreases amplitude more. it also requires more computational power so order of 4 is good enough
[b,a] = butter(N,Fc,"low");
ecg_butter_filtered = filter(b,a,ecg_hfn);

noise_content = ecg_hfn - ecg_butter_filtered;
n_l = length(noise_content);
noise_power = sum(abs(noise_content).^2) /n_l;
N = 10; 

ecg_hfns=fftshift(abs(fft(ecg_hfn,length(ecg_hfn))/length(ecg_hfn)));
f=linspace(-Fs/2,Fs/2,length(ecg_butter_filtered));
f=transpose(f);

ecg_hfnd = zeros(length(ecg_hfns),length(ecg_hfns));
ecg_hfnd(1:length(ecg_hfn)) = ecg_hfn;
ecg_denoised_wiener = wiener2(ecg_hfnd, [1, N],noise_power);
ecg_denoised_wiener = ecg_denoised_wiener(1:length(ecg_hfn),1);
ecg_denoised_wiener_s = fftshift(abs(fft(ecg_denoised_wiener,length(ecg_denoised_wiener))/length(ecg_denoised_wiener)));

% Plot the original and denoised signals for comparison
figure;
set(gcf,'color','w');
plot(ecg_hfn);
title('Original Noisy ECG Signal');
xlabel('Samples');
ylabel('Amplitude');
hold on
plot(ecg_denoised_wiener);
title('Denoised ECG Signal using Wiener Filter');
xlabel('Samples');
ylabel('Amplitude');
legend("ecg_hfn","filtered ecg_hfn")

figure;
set(gcf,'color','w');
plot(f,ecg_hfns)
hold on
plot(f,ecg_denoised_wiener_s)    
title("Frequency spectrum of ecg signal and its filtered version with Wiener")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_hfn","filtered ecg_hfn")

%% Adaptive Filter With LMS (Least Mean Square) Algorithm 
clc;
clear all;
close all;
load('ecg_1.mat');
load('ecg_2.mat');
load('ecg_3.mat');

N = length(ecg_hfn);
Fs = 1000; % sample frequency

% Parameters for LMS algorithm
mu = 0.0001; % Step size (learning rate)
M = 10;    % Order of the filter
w = zeros(M, 1); % Initial filter weights

% Initialize variables for storing the filtered signal
ecg_filtered = zeros(N, 1);

% Adaptive LMS filtering
for n = M:N
    % Input vector from the signal
    x = ecg_hfn(n:-1:n-M+1);
  
    % Output of the filter
    y = w' * x;
    
    % Error calculation (if the reference signal is available)
    % e = reference_signal(n) - y; 
    % For noise cancellation, this is usually the noisy signal itself
    e = ecg_hfn(n) - y; 
    
    % Update filter weights
    w = w + mu * e * x;
    
    % Store the output
    ecg_filtered(n) = y;
end

% Plot the original and filtered signals for comparison
figure;
set(gcf,'color','w');
plot(ecg_hfn);
title('Original Noisy ECG Signal');
xlabel('Samples');
ylabel('Amplitude');
hold on
plot(ecg_filtered);
title('Denoised ECG Signal using Adaptive LMS Filter');
xlabel('Samples');
ylabel('Amplitude');
legend("ecg_hfn","filtered ecg_hfn")

ecg_hfn_freq = abs(fftshift(fft(ecg_hfn,N)))/N;
ecg_hfn_filtered_freq = abs(fftshift(fft(ecg_filtered,N)))/N;
f = linspace(-Fs/2,Fs/2,N);

figure;
set(gcf,'color','w');
plot(f, ecg_hfn_freq);
hold on;
plot(f, ecg_hfn_filtered_freq);
grid on;
xlabel("frequency (Hz)");
ylabel("magnitude");
title('Original and Filtered ECG Signals in frequency domain');
legend("original ecg signal","filtered ecg_lfn signal");
