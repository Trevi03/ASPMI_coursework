clear all
close all
clc;

addpath Algorithms/
%%
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
L = 1024;
w = (0:(L-1)) .* (fs / L);

% x is now the collection of basis vectors defining the DFT
x = (1/L)*exp(1j*2*(0:N-1)'*pi*(0:(L-1))/L).';

mu = 1;
gamma = [0 0.01 0.05 0.5]; % leakage

for i = 1:length(gamma)
    % implement clms algorithm 
    [a,error] = DFT_CLMS(x,y,mu,L,gamma(i));

    H = abs(a);

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    % Plot time-frequency diagram
    subplot(2,2,i)
    surf(1:N,w, H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    ylim([0 700])
    title(sprintf('$\\gamma$ = %0.3f',gamma(i)),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    hold on
end

sgtitle('Time Frequency Estimation using DFT-CLMS','Interpreter','latex','fontsize',20)


%% 
clear
load("EEG_Data/EEG_Data_Assignment1.mat")

N = 1200;
n = 1:N;
v = 0.05; % noise variance

% generate y(n)
a = 1000;
y = POz(a:a+N-1);
y =  y-mean(y);

L = 1024;
w = (0:(L-1)) .* (fs / L);

% x is now the collection of basis vectors defining the DFT
x = (1/L)*exp(1j*2*(0:N-1)'*pi*(0:(L-1))/L).';

colour = {'r','b','c','m','g','y'};

mu = 1;
gamma = [0 0.01 0.05 0.5]; % leakage

for i = 1:length(gamma)
    % implement clms algorithm 
    [a,error] = DFT_CLMS(x,y',mu,L,gamma(i));

    H = abs(a);

    % Remove outliers in the matrix H
    medianH = 50*median(median(H));
    H(H > medianH) = medianH;

    % Plot time-frequency diagram
    subplot(2,2,i)
    surf(1:N,w, H, 'LineStyle','none');
    view(2)
    colorbar('TickLabelInterpreter', 'latex')
    xlabel('Time (n)','fontsize',15)
    ylabel('Frequency (Hz)','fontsize',15)
    ylim([0 100])
    title(sprintf('$\\gamma$ = %0.3f',gamma(i)),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    grid on
    grid minor
    hold on
end

sgtitle('ECG Time Frequency Estimation using DFT-CLMS','Interpreter','latex','fontsize',20)


