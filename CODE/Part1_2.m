clear all
close all
clc;

%% Part (a)
load sunspot.dat

sunData = sunspot(:,2);

% remove mean and detrend
sunDataMean = sunData - mean(sunData);
sunDataDT = detrend(sunDataMean);

% log and subtract mean
sunDataLog = log(sunData + eps) - mean(log(sunData + eps));

% periodogram using a rectangular window
dataLen = ones(1,length(sunDataMean));
[psdOrig,wOrig] = periodogram(sunData,dataLen,[],1);
[psdMean,wMean] = periodogram(sunDataMean,dataLen,[],1);
[psdDT,wDT] = periodogram(sunDataDT,dataLen,[],1);
[psdLog,wLog] = periodogram(sunDataLog,dataLen,[],1);

psdOrig(psdOrig<0.001) = 1;
psdMean(psdMean<0.001) = 1;
psdDT(psdDT<0.001) = 1;
psdLog(psdLog<0.001) = 1;

% Plot
figure;
subplot(1,2,1);
plot(wOrig,10*log10(psdOrig),'b', 'Linewidth', 1.5)
hold on
plot(wMean,10*log10(psdMean),'r', 'Linewidth', 1.5)
hold on
plot(wDT,10*log10(psdDT),'m',  'Linewidth', 1.5)
plot(wLog,10*log10(psdLog),'c', 'Linewidth',1.5)
ax = gca;
ax.FontSize = 14;
xlabel('Normalized frequency')
ylabel('PSD (dB) [10*log_{10}(X)]')
title('Periodogram Sunspot Data Series - Rectangular Window', 'Fontsize', 16)
legend('Raw Data','Mean Removal','Detrend','Log + Mean Removal','fontsize',13.5)
grid on
grid minor


% periodogram using a hamming window
[psdOrig,wOrig] = periodogram(sunData,dataLen,[],1);
[psdMean,wMean] = periodogram(sunDataMean,hamming(length(sunDataMean)),[],1);
[psdDT,wDT] = periodogram(sunDataDT,hamming(length(sunDataDT)),[],1);
[psdLog,wLog] = periodogram(sunDataLog,hamming(length(sunDataLog)),[],1);

psdOrig(psdOrig<0.001) = 1;
psdMean(psdMean<0.001) = 1;
psdDT(psdDT<0.001) = 1;
psdLog(psdLog<0.001) = 1;

% Plot
subplot(1,2,2);
plot(wOrig,10*log10(psdOrig),'b', 'Linewidth', 1.5)
hold on
plot(wMean,10*log10(psdMean),'r', 'Linewidth', 1.5)
hold on
plot(wDT,10*log10(psdDT),'m',  'Linewidth', 1.5)
plot(wLog,10*log10(psdLog),'c', 'Linewidth',1.5)
ax = gca;
ax.FontSize = 14;
xlabel('Normalized frequency')
ylabel('PSD (dB) [10*log_{10}(X)]')
title('Periodogram Sunspot Data Series - Hamming Window', 'Fontsize', 16)
legend('Raw Data','Mean Removal','Detrend','Log + Mean Removal','fontsize',13.5)
grid on
grid minor


%% Part (b)
load ../EEG_Data/EEG_Data_Assignment1.mat

N = length(POz);
winSz = [1,5,10]*fs;

% remove mean
POz = POz - mean(POz);

% Standard PSD
[psdStan,fStan] = pwelch(POz, ones(1,N), 0,fs*10, fs, 'onesided');

% Welch's method - Averaged overlapping windows
[psd1, f1] = pwelch(POz, rectwin(winSz(1)),0,fs*10, fs, 'onesided');
[psd5, f5] = pwelch(POz, rectwin(winSz(2)),0,fs*10, fs, 'onesided');
[psd10, f10] = pwelch(POz, rectwin(winSz(3)),0,fs*10, fs, 'onesided');  


% Plot
figure;
subplot(2,2,1);
plot(fStan, 10*log10(psdStan),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('Standard EEG Periodogram','fontsize',15)
grid on
grid minor

subplot(2,2,2)
plot(f1,10*log10(psd1),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 1s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,3)
plot(f5,10*log10(psd5), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 5s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,4)
plot(f10,10*log10(psd10), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 10s Windowing','fontsize',15)
grid on
grid minor

sgtitle('Averaged Periodogram for EEG Data - Rectangular Window', 'Fontsize', 20)
set(gcf,'color','w')



% Welch's method - Hamming
[psd1, f1] = pwelch(POz, hamming(winSz(1)),0,fs*10, fs, 'onesided');
[psd5, f5] = pwelch(POz, hamming(winSz(2)),0,fs*10, fs, 'onesided');
[psd10, f10] = pwelch(POz, hamming(winSz(3)),0,fs*10, fs, 'onesided');  

% Plot
figure;
subplot(2,2,1);
plot(fStan, 10*log10(psdStan),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('Standard EEG Periodogram','fontsize',15)
grid on
grid minor

subplot(2,2,2)
plot(f1,10*log10(psd1),'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 1s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,3)
plot(f5,10*log10(psd5), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD (dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 5s Windowing','fontsize',15)
grid on
grid minor

subplot(2,2,4)
plot(f10,10*log10(psd10), 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('EEG Periodogram - 10s Windowing','fontsize',15)
grid on
grid minor

sgtitle('Averaged Periodogram for EEG Data - Hamming Window', 'Fontsize', 20)
set(gcf,'color','w')

%% Just to illustrate variance and bias
figure;
plot(fStan, 10*log10(psdStan),'Color','blue','LineWidth',1,DisplayName="Standard")
hold on
plot(f10,10*log10(psd1), 'Color','red','LineWidth',3,DisplayName="1s Windowing")
ax = gca;
ax.FontSize = 15;
xlabel('Normalized frequency')
ylabel('PSD(dB)')
xlim([0,100])
ylim([-150, -80])
title('Standard vs Averaged Periodogram','fontsize',15)
legend
grid on
grid minor
set(gcf,'color','w')
