clc, clear, close all;
load('ecg_1.mat')
load('ecg_2.mat')
load('ecg_3.mat')

%butterworth 
Fs=1000;
Ts=1/Fs;
Fc=50/Fs; %cutoff freq/Fs for normalization. cutoff freq determined by observing freq spectrum of the ecg signal
N=4; %filter order
% by increasing order we get a sharper filter freq response supresses freqs above Fc better. 
% but increases phase shift and decreases amplitude more. it also requires more computational power so order of 4 is good enough
[b,a]=butter(N,Fc,"low");
ecg_butter_filtered=filter(b,a,ecg_hfn);

figure (1)
plot(ecg_hfn)
hold on
plot(ecg_butter_filtered)
title("ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_2","filtered ecg_2")

ecg_hfns=fftshift(abs(fft(ecg_hfn,length(ecg_hfn))/length(ecg_hfn)));
ecg_butter_filtereds=fftshift(abs(fft(ecg_butter_filtered,length(ecg_butter_filtered))/length(ecg_butter_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_butter_filtered));
f=transpose(f);

figure (2)
plot(f,ecg_hfns)
hold on
plot(f,ecg_butter_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

waveletType = 'db4'; %sym4
noise_content=ecg_hfn-ecg_butter_filtered;
variance=var(noise_content);

%multilevel
level=5;
[c,l]=wavedec(ecg_hfn,level,waveletType);
denoised_ecg_wavelet=[];

j=1;
len=0;
for i=level:-1:1
    cD = detcoef(c,l, i); %detail coef in each level
    th=sqrt(2*variance*log(length(cD)))   
    cD = wthresh(cD, 's', th);
    len=len+l(j);
    c(len+1:len+length(cD))=cD;
    j=j+1;
end

denoised_ecg_wavelet=waverec(c,l,waveletType);
denoised_ecg_wavelet_s=fftshift(abs(fft(denoised_ecg_wavelet,length(denoised_ecg_wavelet))/length(denoised_ecg_wavelet)));


% Plot original and denoised ECG signal

figure(3)
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
figure (4)
plot(f,ecg_hfns)
hold on
plot(f,denoised_ecg_wavelet_s)    
title("frequency spectrum of ecg signal and its filtered version with wavelet")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

