clear all
close all
clc;

addpath Algorithms/

%%

N = 1500;
fs = 2000;
v = 0.05; % noise variance

% circular complex-valued white noise
eta = sqrt(v).*randn(1,N) + 1j*sqrt(v).*randn(1,N); 
% dphi/dn
phi_dot = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
phi = cumtrapz(phi_dot);

% plot the f(n) and phi(n) for demonstration
figure
subplot(1,2,1)
plot(phi_dot,'LineWidth',2,Color="#77AC30")
xlabel('Time (n)')
ylabel('Frequency (Hz)')
title('Frequency f(n)','Interpreter','latex')
ax = gca;
ax.FontSize = 15;
grid on
grid minor

subplot(1,2,2)
plot(wrapTo2Pi(phi),'LineWidth',2,Color="#A2142F")
xlabel('Time (n)')
xlim([0 200])
ylabel('Angle (rads)')
title('Phase $\Phi(n)$','Interpreter','latex')
ax = gca;
ax.FontSize = 15;
grid on
grid minor
set(gcf,'color','w')

%%
y = exp(1j*((2*pi)/fs)*phi) + eta;

orders = [1,5,10];
colour = {'r','b','c','m','g','y'};

figure
hold on
for i = 1:length(orders)
    % using aryule to find the AR(1) coefficient
    a = aryule(y, orders(i));

    % N-point frequency response vector h, and angular frequency vector w 
    [h,w] = freqz(1,a,N,fs);

    % compute the psd
    psd = 10*log10(abs(h).^2);
    subplot(1,3,i)
    plot(w,psd,'color',colour{i},'LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)','fontsize',10)
    title(sprintf('AR(%0.0f)',orders(i)),'Interpreter','latex','fontsize',18)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
end
sgtitle('Power Spectra of signal y(n)','fontsize',20)

% splitting into 500 segment lengths in order to plot the frequency
figure
hold on
for i = 1:3
    % using aryule to find the AR(1) coefficinet for the complete signal
    a = aryule(y(500*(i-1)+1:500*i), 1);
    % obtain the N-point frequency response vector, and the corresponding
    % angular frequency vector w - for the digitial filter constructed
    [h,w] = freqz(1,a,N/3,fs);
    % compute the power spectral density
   
    psd = 10*log10(abs(h).^2);
    subplot(1,3,i)
    plot(w,psd,colour{i},'LineWidth',2)
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)','fontsize',10)
    title(sprintf('AR(1) Segment %0.0f',i),'Interpreter','latex','fontsize',18)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    
end
sgtitle('Power Spectra of y(n) with 3 segments of length 500','fontsize',20)

%% CLMS
clear

N = 1500;
n = 1:N;
fs = 2000;
v = 0.05; % noise variance

% circular complex-valued white noise
eta = sqrt(v).*randn(1,N) + 1j*sqrt(v).*randn(1,N); 
% dphi/dn
phi_dot = [100*ones(1,500), 100 + ((501:1000)-500)/2, 100 + (((1001:1500)-1000)/25).^2];
phi = cumtrapz(phi_dot);
% generate y(n)
y = exp(1j*((2*pi)/fs)*phi) + eta;

colour = {'r','b','c','m','g','y'};

mu = [0.001 0.01 0.1 0.5];

figure
hold on
% shift input array in preparation for prediction
x = [0, y(1:end-1)];
for i = 1:length(mu)
    % pre allocation for coefficinets, error and signal
    a = complex(zeros(1,N));  % need to change notation for psd code part
    error = complex(zeros(1,N));

    % implement clms algorithm 
    [a,error] = predCLMS(x,y,mu(i),1);

    L = 1024;
    H = zeros(L, N);

    % code provided in the coursework documentation
    for n=1:N
        % Run complex-valued LMS algorthm to esitmate AR coefficient a_hat1(n)
        [h, w] = freqz(1, [1; -conj(a(n))], L);
        H(:, n) = abs(h).^2;
    end

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    % Plot time-frequency diagram
    subplot(2,2,i)
    surf(1:N,((w*fs)/(2*pi)), H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time Index (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    title(sprintf('$\\mu$ = %0.3f',mu(i)),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    hold on
end

sgtitle('Time Frequency Estimation using CLMS','Interpreter','latex','fontsize',20)





