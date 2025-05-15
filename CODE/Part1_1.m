clear all
close all
clc;

%%

N = 2e3;        % number of discrete samples
f = 1e3;        % frequency
T = 1/f;        % period
time = 1:N;     % time vector array

% Impulse
impulse = zeros(1,N);
impulse(N/2) = 30;

% Fast ACF decay
sinusoid = cos(2*pi*T*time);

% Slow ACF decay
alpha = 0.01;
dampedSinusoid = exp(-alpha * time) .* cos(2*pi*T*time);


% ACF
acfImpulse = xcorr(impulse,'biased'); 
acfSine = xcorr(sinusoid,'biased');     % sine wave acf
acfDamp = xcorr(dampedSinusoid,'biased');  % unit impulse acf

% PSD def 1
psdImpulse1 = abs(fftshift(fft(acfImpulse)));
psdDamp1 = abs(fftshift(fft(acfSine)));
psdSine1 = abs(fftshift(fft(acfDamp)));

% PSD def 2
psdImpulse2 = abs(fftshift(fft(impulse))).^2 / N;
psdDamp2 = abs(fftshift(fft(dampedSinusoid))).^2 / N;
psdSine2 = abs(fftshift(fft(sinusoid))).^2 / N;

% frequency axis in Hz for psd
freqRad1 = [-pi+pi/length(acfSine):2*pi/length(acfSine):pi-pi/length(acfSine)];
freqHz1 = freqRad1./(2*pi).*f;

freqRad2 = [-pi+pi/length(psdSine2):2*pi/length(psdSine2):pi-pi/length(psdSine2)];
freqHz2 = freqRad2./(2*pi).*f;


figure;
subplot(3, 1, 1);
plot([-N+1:N-1],acfSine,'b','LineWidth', 2);
hold on;
plot([-N+1:N-1],acfImpulse,'r', 'LineWidth', 2);
legend('Sinusoid','Impulse');
ax = gca;
ax.FontSize = 15;
title('Signal ACFs','fontsize',15)
xlabel('Lags (sample)');
ylabel('ACF');

grid on
grid minor


subplot(3, 1, 2);
plot(freqHz1, psdSine1,'b', 'LineWidth', 2);
hold on;
plot(freqHz2, psdSine2,'r--', 'LineWidth', 2);
legend('Def 1', 'Def 2');
title('Periodogram of Sinusoid','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Normalised frequency (\pi rads/sample)');
ylabel('PSD');
grid on 
grid minor


subplot(3, 1, 3);
plot(freqHz1, psdImpulse1,'b', 'LineWidth', 2);
hold on;
plot(freqHz2, psdImpulse2, 'r--', 'LineWidth', 2);
legend('Def 1', 'Def 2');
title('Periodogram of Impulse','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Normalised frequency (\pi rads/sample)');
ylabel('PSD');
grid on
grid minor
set(gcf,'color','w')

