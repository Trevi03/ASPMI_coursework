clear all
close all
clc;

%% Part (a)
N = 1000;
time = 0:N-1;
WGN = randn(1,N);
nSine = sin(time*0.1) + WGN;
filtWGN = filter([1/8 1/8 1/8 1/8 1/8 1/8 1/8 1/8], 1, WGN); % Moving Average Filter

% calculate biased  ACF
[WGN_b, WGN_blag] = xcorr(WGN, 'biased');
[nSine_b, nSine_blag] = xcorr(nSine, 'biased');
[filtWGN_b, filtWGN_blag] = xcorr(filtWGN, 'biased');

% calculate unbiased  ACF
[WGN_u, WGN_ulag] = xcorr(WGN, 'unbiased');
[nSine_u, nSine_ulag] = xcorr(nSine, 'unbiased');
[filtWGN_u, filtWGN_ulag] = xcorr(filtWGN, 'unbiased');

% make 0 centered
psdWGN_b = fftshift(real(fft(ifftshift(WGN_b))));
psdNSine_b = fftshift(real(fft(ifftshift(nSine_b))));
psdFilt_b = fftshift(real(fft(ifftshift(filtWGN_b))));

psdWGN_u = fftshift(real(fft(ifftshift(WGN_u))));
psdNSine_u = fftshift(real(fft(ifftshift(nSine_u))));
psdFilt_u = fftshift(real(fft(ifftshift(filtWGN_u))));


%%

% Plot WGN
figure;
subplot(3,2,1)
plot(WGN_ulag, WGN_u, 'Color', 'blue','LineWidth',2)
hold on
plot(WGN_blag, WGN_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('ACF')
title('WGN ACF')
legend('Unbiased','Biased')
grid on
grid minor

subplot(3,2,2)
plot(WGN_ulag, psdWGN_u, 'Color', 'blue','LineWidth',2)
hold on
plot(WGN_blag, psdWGN_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('PSD')
title('WGN Correlogram')
legend('Unbiased','Biased')
grid on
grid minor

% noisy Sinusoid
subplot(3,2,3)
plot(nSine_ulag, nSine_u, 'Color', 'blue','LineWidth',2)
hold on
plot(nSine_blag, nSine_b, 'Color','red','LineWidth',2) 
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('ACF')
title('Noisy Sinusoid ACF')
legend('Unbiased','Biased')
grid on
grid minor


subplot(3,2,4)
plot(nSine_ulag, psdNSine_u, 'Color', 'blue','LineWidth',2)
hold on
plot(nSine_blag, psdNSine_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('PSD')
title('Noisy Sinusoid Correlogram')
legend('Unbiased','Biased')
grid on
grid minor

% Filtered WGN
subplot(3,2,5)
plot(filtWGN_ulag, filtWGN_u, 'Color', 'blue','LineWidth',2)
hold on
plot(filtWGN_blag, filtWGN_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('ACF')
title('Filtered WGN ACF')
legend('Unbiased','Biased')
grid on
grid minor

subplot(3,2,6)
plot(filtWGN_ulag, psdFilt_u, 'Color', 'blue','LineWidth',2)
hold on
plot(filtWGN_blag, psdFilt_b, 'Color','red','LineWidth',2)
ax = gca;
ax.FontSize = 13;
xlabel('Lag Values')
xlim([-N, N])
ylabel('PSD')
title('Filtered WGN Correlogram')
legend('Unbiased','Biased')
grid on
grid minor

%% Part (b)

fs = 10;
N = 1000;
n = 0:1/(2*fs):fs;
x = [1.5*sin(2*pi*n) + 0.5*sin(2*pi*n*0.3) + 1.2*sin(2*pi*n*1.8) zeros(1, N-length(n))];
all_psds = zeros(100, N*2-1);


figure;
hold on
for i = 1:100
    noisy_signal = x + randn(1,length(x)); %wgn(length(test_signal), 1, 1)';
    % Correlogram
    [corr_b, corr_lag] = xcorr(noisy_signal, 'biased');
    all_psds(i,:) = fftshift(real(fft(ifftshift(corr_b))));
    subplot(1,2,1)
    f_rad = (corr_lag/max(corr_lag))*fs;
    plot(f_rad, real(all_psds(i,:)),'color','c','LineWidth',1.5);
    hold on
end

f_rad = (corr_lag/max(corr_lag))*fs;
% Mean
plot(f_rad, mean(real(all_psds)),'color','b','LineWidth',2);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency [Hz]')
ylabel('PSD')
set(gca,'fontsize',15)
title('PSD estimates (different realisations and mean)','fontsize',16)
grid on 
grid minor

subplot(1,2,2)
% Standard Deviation
plot(f_rad, std(real(all_psds)),'color','r','LineWidth',2);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency [Hz]')
ylabel('PSD')
set(gca,'fontsize',15)
title('Standard deviation of the PSD estimate','fontsize',16)
grid on
grid minor
set(gcf,'color','w')


%% Part (c)

figure;
hold on
for i = 1:100
    subplot(1,2,1)
    plot(f_rad, 10*log10(real(all_psds(i,:))),'color','c','LineWidth',1.5);
    hold on
end

% Mean
plot(f_rad, 10*log10(mean(real(all_psds))),'color','b','LineWidth',2);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency [Hz]')
ylabel('PSD (dB)')
set(gca,'fontsize',15)
title('PSD estimates (different realisations and mean)','fontsize',16)
grid on 
grid minor

subplot(1,2,2)
% Standard Deviation
plot(f_rad, std(10*log10(real(all_psds))),'color','r','LineWidth',2);
xlim([0,3])
ax = gca;
ax.FontSize = 15;
xlabel('Frequency [Hz]')
ylabel('PSD (dB)')
set(gca,'fontsize',15)
title('Standard deviation of the PSD estimate','fontsize',16)
grid on
grid minor
set(gcf,'color','w')


%% Part (d)
N= 1000;
len = [30 35 40 45 50];
colors = {'y','r','g','b','m'};

figure;
hold on
for i = 1: length(len)
    n = 0:len(i);
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    zeroArr = zeros(1, N-length(x));
    psd_data = abs((fft([x zeroArr]))./length(n));
    f_axis = [0:N-1]/N;
    plot(f_axis, 10*log10(psd_data),'color',colors{i},'LineWidth',2)
    hold on
end
ax = gca;
ax.FontSize = 14;
xlabel('Frequency (Hz)')
ylabel('PSD (dB)')
title('PSD: Noisy Complex Exponential Signal','fontsize',14)
xlim([0.2 0.5])
legend('N=30','N=35','N=40','N=45','N=50')
grid on
grid minor


%% Part (e)

N = 256;
n = 0:30-1;
all_psds = zeros(100, N);


figure;
subplot(1,2,1);
hold on
for i = 1:100
    noise = 0.2/sqrt(2)*(randn(size(n))+1j*randn(size(n)));
    x = exp(1j*2*pi*0.3*n) + exp(1j*2*pi*0.32*n) + noise;
    [~,R] = corrmtx([x zeros(1, N-length(x))],14,'mod');
    [all_psds(i, :), F] = pmusic(R,2,[ ],1,'corr');
    plot(F,all_psds(i, :),'color','c','linewidth',2); 
end

% Mean
plot(F, mean(all_psds),'color','b','LineWidth',2);
ax = gca;
ax.FontSize = 15;
set(gca,'xlim',[0.25 0.40]);
set(gca,'fontsize',15);
xlabel('Hz'); 
ylabel('Pseudospectrum');
title('MUSIC Pseudospectrum - Realisations and mean','fontsize',16)
grid on 
grid minor

subplot(1,2,2)
% Standard Deviation
plot(F, std(all_psds),'color','r','LineWidth',2);
ax = gca;
ax.FontSize = 15;
xlabel('Frequency [\pi radians]')
ylabel('PSD')
set(gca,'fontsize',15)
set(gca,'xlim',[0.25 0.40]);
title('Standard deviation of MUSIC Pseudospectrum','fontsize',16)
grid on
grid minor
set(gcf,'color','w')
