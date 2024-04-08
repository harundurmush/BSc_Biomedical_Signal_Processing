clc, clear, close all;
%% QUESTION 1 
%a
load('ecg_1.mat')
load('ecg_2.mat')
load('ecg_3.mat')
Fs=1000;
Ts=1/Fs;
ecg_filtered=zeros(length(ecg_hfn),1);
plot(ecg_hfn);
for i=1:length(ecg_hfn)
    if i<5
        ecg_filtered(i)=mean(ecg_hfn(1:i+5)); %need special cases for the first and last indexes of signal array
    elseif i>length(ecg_hfn)-5
        ecg_filtered(i)=mean(ecg_hfn(i:length(ecg_hfn)));
    else
        ecg_filtered(i)=mean(ecg_hfn(i-4:i+5)); % to achieve a window size of 10
    end
end

figure (1)
plot(ecg_hfn)
hold on
plot(ecg_filtered)
title("ecg signal and its filtered version with moving average filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_2","filtered ecg_2")

%b
ecg_hfns=fftshift(abs(fft(ecg_hfn,length(ecg_hfn))/length(ecg_hfn)));
ecg_filtereds=fftshift(abs(fft(ecg_filtered,length(ecg_filtered))/length(ecg_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_filtered));
f=transpose(f);
figure (2)
plot(f,ecg_hfns)
hold on
plot(f,ecg_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with moving average filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

%% QUESTION 2
%butterworth
%better filtering than moving avg filter but more loss on amplitude
Fc=50/Fs; %cutoff freq/Fs for normalization. cutoff freq determined by observing freq spectrum of the ecg signal
N=4; %filter order
% by increasing order we get a sharper filter freq response supresses freqs above Fc better. 
% but increases phase shift and decreases amplitude more. it also requires more computational power so order of 4 is good enough
[b,a]=butter(N,Fc,"low");
ecg_butter_filtered=filter(b,a,ecg_hfn);
figure (3)
plot(ecg_hfn)
hold on
plot(ecg_butter_filtered)
title("ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_2","filtered ecg_2")

ecg_butter_filtereds=fftshift(abs(fft(ecg_butter_filtered,length(ecg_butter_filtered))/length(ecg_butter_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_butter_filtered));
f=transpose(f);
figure (4)
plot(f,ecg_hfns)
hold on
plot(f,ecg_butter_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with 4th order butterworth filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

%% Chebychev type1 
N=4;
Rp=1; %ripple dB in passband
[b,a]=cheby1(N,Rp,Fc,"low");

ecg_cheby1_filtered=filter(b,a,ecg_hfn);
amplitude=max(ecg_cheby1_filtered)

figure (5)
plot(ecg_hfn)
hold on
plot(ecg_cheby1_filtered)
title("ecg signal and its filtered version with 4th order chebyschev type1 filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_2","filtered ecg_2")

ecg_cheby1_filtereds=fftshift(abs(fft(ecg_cheby1_filtered,length(ecg_cheby1_filtered))/length(ecg_cheby1_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_cheby1_filtered));
f=transpose(f);
figure (6)
plot(f,ecg_hfns)
hold on
plot(f,ecg_cheby1_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with 4th order chebyschev type1 filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

%% Chebychev type2
N=8; %higher order reduced amplitude loss in contrast to other filter designs
Rs=20; %ripple dB in passband
%as i increase the dB it provides a better filtering but the loss in amplitude is much higher than cheby1 and butter
[b,a]=cheby2(N,Rs,Fc,"low");

ecg_cheby2_filtered=filter(b,a,ecg_hfn);
figure (7)
plot(ecg_hfn)
hold on
plot(ecg_cheby2_filtered)
title("ecg signal and its filtered version with 8th order chebyschev type2 filter")
xlabel("Sample");
ylabel("Amplitude");
legend("ecg_2","filtered ecg_2")

ecg_cheby2_filtereds=fftshift(abs(fft(ecg_cheby2_filtered,length(ecg_cheby2_filtered))/length(ecg_cheby2_filtered)));
f=linspace(-Fs/2,Fs/2,length(ecg_cheby2_filtered));
f=transpose(f);
figure (8)
plot(f,ecg_hfns)
hold on
plot(f,ecg_cheby2_filtereds)   
title("frequency spectrum of ecg signal and its filtered version with 8th order chebyschev type2 filter")
xlabel("Frequency(Hz)");
ylabel("Magnitude");
legend("ecg_2","filtered ecg_2")

%% QUESTION 3
figure (9)
stft(ecg_hfn,Fs,"Window",rectwin(256),"OverlapLength",128);
title("STFT with rectangular window")
figure (10)
stft(ecg_hfn,Fs,"Window",hamming(256),"OverlapLength",128);
title("STFT with hamming window")

%% QUESTION 4
%periodogram method that uses rectangular window
[Pxx,f]=periodogram(ecg_lfn);
figure (11)
subplot(3,1,1)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_1 with periodogram")
[Pxx,f]=periodogram(ecg_hfn);
subplot(3,1,2)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_2 with periodogram")
[Pxx,f]=periodogram(ecg_noisy);
subplot(3,1,3)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_3 with periodogram")


%% modified periodogram with hamming window
[Pxx,f]=periodogram(ecg_lfn,hamming(length(ecg_lfn)),length(ecg_lfn),Fs);
figure (14)
subplot(3,1,1)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_1 with periodogram with hamming window")
[Pxx,f]=periodogram(ecg_hfn,hamming(length(ecg_hfn)),length(ecg_hfn),Fs);
subplot(3,1,2)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_2 with periodogram with hamming window")
[Pxx,f]=periodogram(ecg_noisy,hamming(length(ecg_noisy)),length(ecg_noisy),Fs);
subplot(3,1,3)
plot(f,10*log10(Pxx)) %log for dB scale
title("PSD of ecg_3 with periodogram with hamming window")
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");

%% welch method
[Pxx,f]=pwelch(ecg_lfn,rectwin(length(ecg_lfn)),[],length(ecg_lfn),Fs);
figure (17)
subplot(3,1,1)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_1 with periodogram with welch method")
[Pxx,f]=pwelch(ecg_hfn,rectwin(length(ecg_hfn)),[],length(ecg_hfn),Fs);
subplot(3,1,2)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_2 with periodogram with welch method")
[Pxx,f]=pwelch(ecg_noisy,rectwin(length(ecg_noisy)),[],length(ecg_noisy),Fs);
subplot(3,1,3)
plot(f,10*log10(Pxx)) %log for dB scale
ylabel("Power/Frequency(dB/Hz)");
xlabel("Frequency(Hz)");
title("PSD of ecg_3 with periodogram with welch method")

%% levinson-durbin method
[r,lags] = xcorr(ecg_hfn, 'biased');
[a,e]=levinson(r,10);
[h,f]=freqz(sqrt(e),a,length(ecg_hfns),Fs);
psd=10*log10(abs(h));
figure (20)
plot(f,psd);
xlabel('Frequency(Hz)');
ylabel('Power/Frequency(dB/Hz)');
title('PSD of ecg_2 using levinson-durbin algorithm');
ecg_lfns=fftshift(abs(fft(ecg_lfn,length(ecg_lfn))/length(ecg_lfn)));
[r,lags] = xcorr(ecg_lfn, 'biased');
[a,e]=levinson(r,10);
[h,f]=freqz(sqrt(e),a,length(ecg_lfns),Fs);
psd=10*log10(abs(h));
figure (21)
plot(f,psd);
xlabel('Frequency(Hz)');
ylabel('Power/Frequency(dB/Hz)');
title('PSD of ecg_2 using levinson-durbin algorithm');
ecg_noisys=fftshift(abs(fft(ecg_noisy,length(ecg_noisy))/length(ecg_noisy)));
[r,lags] = xcorr(ecg_noisy, 'biased');
[a,e]=levinson(r,10);
[h,f]=freqz(sqrt(e),a,length(ecg_noisys),Fs);
psd=10*log10(abs(h));
figure (22)
plot(f,psd);
xlabel('Frequency(Hz)');
ylabel('Power/Frequency(dB/Hz)');
title('PSD of ecg_2 using levinson-durbin algorithm');

%% QUESTION 5
%coefficients are found from the equations in the lecture slide
%lowpass filter
b = [1,0,0,0,0,0,-2,0,0,0,0,0,1]; 
a = [1,-2,1]; 
ecg_lpf=filter(b,a,ecg_hfn);
%highpass filter
b = [-1/32,zeros(1, 15),1,-1,zeros(1, 15),1/32]; 
a = [1,-1];
ecg_hpf=filter(b,a,ecg_lpf);
b=[2,1,0,-1,-2]; 
a=8;
ecg_dif = filter(b,a,ecg_hpf);
%squaring
for i=1:length(ecg_dif)
    ecg_dif(i)=ecg_dif(i)^2;
end
%moving average filter 
%using same design in Q1
ecg_qrs=zeros(length(ecg_dif),1);
for i=1:length(ecg_qrs)
    if i<5
        ecg_qrs(i)=mean(ecg_dif(1:i+5)); 
    elseif i>length(ecg_dif)-5
        ecg_qrs(i)=mean(ecg_dif(i:length(ecg_dif)));
    else
        ecg_qrs(i)=mean(ecg_dif(i-4:i+5)); 
    end
end
ecg_qrs=ecg_qrs./20; %sclaed the signal down for easier obeservation on the plot. amplitude was too high for comparison
figure (23)
plot(ecg_hfn)
hold on
plot(ecg_qrs)
legend("ecg_2","QRS")
title("ecg signal and QRS part detected signal")
xlabel("Sample");
ylabel("Amplitude")

