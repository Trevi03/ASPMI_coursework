clear all
close all
clc;

addpath Algorithms/

%% Part (a)
N = 1000;
nTrials = 100;
deltas = [1,2,3,4];
mu = 0.01; 
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2) (+0v(n-1))

% generating the clean s(n)
x = sin(0.01*pi.*(0:N-1)); 

figure
hold on
for i = 1: length(deltas)
    xhat_all = zeros(nTrials,N);
    s_all =  zeros(nTrials,N);
    MSPE = zeros(1,nTrials); % Mean Square Prediction Error
    
    for j = 1: nTrials
        v = randn(1,N); % unit variance

        % generating the coloured noise component of s(n)
        eta = filter(b,a,v);

        % constructing the input signal to ALE LMS
        s = x + eta;
        s_all(j,:) = s;
        
        [xhat,w,err] = LMS_ALE(s,mu,M,0,deltas(i)); % (x_input,mu,order,leakage,delay)
        xhat_all(j,:) = xhat;
        MSPE(j) = mean((x(100:end)-xhat(100:end)).^2);
        
    end
    subplot(1,4,i)
    p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
    hold on
    p2 = plot(xhat_all','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
    hold on
    p3 = plot(x,'y','LineWidth',2,'DisplayName','$x(n)$')
    MSPE4title = mean(MSPE);
    title(sprintf('MSPE = %0.3f, $\\Delta$ = %0.0f',MSPE4title,deltas(i)),'Interpreter','latex','fontsize',15)
    ax = gca;
    ax.FontSize = 15;
    xlabel('Sample (n)','fontsize',15)
    ylabel('(AU)','fontsize',15)
    set(groot,'defaultLegendInterpreter','latex');
    allps = [p1(1),p2(1),p3];
    legend(allps)
    grid on
    grid minor
    
end
sgtitle('Empirical Justification of Minimum Delay','fontsize',18)


%% Part (b)

% MSPE against Delay
M = [5,10,15,20];
deltas = 1:25;

MSPE_od =  zeros(length(M),length(deltas));
colors = {'b','r','c','m'};

figure;
subplot(1,2,1)
hold on
for i = 1: length(M)
    for j = 1: length(deltas)
        xhat_all = zeros(nTrials,N);
        s_all =  zeros(nTrials,N);
        MSPE = zeros(1,nTrials); % Mean Square Prediction Error

        for k=1:nTrials
            v = randn(1,N); % unit variance
    
            % generating the coloured noise component of s(n)
            eta = filter(b,a,v);
    
            % constructing the input signal to ALE LMS
            s = x + eta;
            s_all(j,:) = s;
            
            [xhat,w,err] = LMS_ALE(s,mu,M(i),0,deltas(j)); % (x_input,mu,order,leakage,delay)
            xhat_all(j,:) = xhat;
            MSPE(k) = mean((x(100:end)-xhat(100:end)).^2);
        end
        MSPE_od(i,j) = mean(MSPE);        
    end
    plot(MSPE_od(i,:),'Color',colors{i},'LineWidth',2)
    hold on
end
ax = gca;
ax.FontSize = 15;
xlabel('Delay (\Delta)')
ylabel('MSPE')
title('Relationship between Delay and MSPE')
legend('M=5','M=10','M=15','M=20')
grid on
grid minor
set(gcf,'color','w')

% Delay against MSPE 
M = 1:20;
delta = 3;

MSPE_delay =  zeros(1,length(M));
colors = {'b','r','c','m'};

subplot(1,2,2)
hold on
for i = 1: length(M)
    xhat_all = zeros(nTrials,N);
    s_all =  zeros(nTrials,N);
    MSPE = zeros(1,nTrials); % Mean Square Prediction Error

    for k=1:nTrials
        v = randn(1,N); % unit variance

        % generating the coloured noise component of s(n)
        eta = filter(b,a,v);

        % constructing the input signal to ALE LMS
        s = x + eta;
        s_all(k,:) = s;
        
        [xhat,w,err] = LMS_ALE(s,mu,M(i),0,delta); % (x_input,mu,order,leakage,delay)
        xhat_all(k,:) = xhat;
        MSPE(k) = mean((x(100:end)-xhat(100:end)).^2);
    end
    MSPE_delay(i) = mean(MSPE);     
end
plot(MSPE_delay,'Color','r','LineWidth',2)
ax = gca;
ax.FontSize = 15;
xlabel('Filter Order')
ylabel('MSPE')
title('Relationship between Filter Order and MSPE, \Delta=3')
grid on
grid minor
set(gcf,'color','w')

%% Part (c)
clear;

N = 1000;
nTrials = 200;
delta = 3;
mu = 0.01; 
M = 5; % filter order (minimum one that is tested in part b)
a = 1;
b = [1 0 0.5]; % since we have eta(n) = v(n) + 0.5v(n-2)

% generating the clean s(n)
x = sin(0.01*pi.*(0:N-1)); 

figure
hold on

s_all =  zeros(nTrials,N);

xhat_allALE = zeros(nTrials,N);
MSPE_ALE = zeros(1,nTrials);   % Mean Square Prediction Error

xhat_allANC = zeros(nTrials,N);
MSPE_ANC = zeros(1,nTrials);

for k=1:nTrials
    v = randn(1,N); % unit variance
    % generating the coloured noise component of s(n)
    eta = filter(b,a,v);

    % constructing the input signal and noise
    s = x + eta;
    u = 1.5*eta+0.1;
    s_all(k,:) = s;
    
    [xhatALE,~,~] = LMS_ALE(s,mu,M,0,delta); % (x_input,mu,order,leakage,delay)
    xhat_allALE(k,:) = xhatALE;
    MSPE_ALE(k) = mean((x(100:end)-xhatALE(100:end)).^2);

    [~,xhatANC,~] = LMS_ANC(s,u,mu,M,0); 
    % [noiseEst,w,xhatANC] = ANC_LMS(s,u,mu,0,M);
    xhat_allANC(k,:) = xhatANC;
    MSPE_ANC(k) = mean((x(100:end)-xhatANC(100:end)').^2);
end

% ALE
subplot(1,3,1)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allALE','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'g','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ALE);
title(sprintf('ALE: MSPE = %0.3f, $\\Delta$ = %0.0f',MSPE4title,delta),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on
grid minor

% ANC
subplot(1,3,2)
p1 = plot(s_all','b','LineWidth',1.5,'DisplayName','$s(n)$')
hold on
p2 = plot(xhat_allANC','r','LineWidth',1.5,'DisplayName','$\hat{x}(n)$')
hold on
p3 = plot(x,'g','LineWidth',2,'DisplayName','$x(n)$')
MSPE4title = mean(MSPE_ANC);
title(sprintf('ANC: MSPE = %0.3f, $\\Delta$ = %0.0f',MSPE4title,0),'Interpreter','latex','fontsize',15)
ax = gca;
ax.FontSize = 15;
xlabel('Sample (n)','fontsize',15)
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
allps = [p1(1),p2(1),p3];
legend(allps)
grid on

% ensemble of realisations
ANC_ensemble = mean(xhat_allANC);
ALE_ensemble = mean(xhat_allALE);
subplot(1,3,3)
plot(ALE_ensemble,'b','LineWidth',2,DisplayName='ALE');
hold on
plot(ANC_ensemble,'r','LineWidth',2,DisplayName='ANC');
hold on
plot(x,'g','LineWidth',2,DisplayName='$x(n)$');
ax = gca;
ax.FontSize = 15;
title('Ensemble Means: ALE vs ANC','Interpreter','latex','fontsize',15)
xlabel('Sample (n)')
ylabel('(AU)','fontsize',15)
set(groot,'defaultLegendInterpreter','latex');
grid on
legend;

%% Part (d)
clear;
load EEG_Data/EEG_Data_Assignment2.mat


% create a time axis
time =  (0:1/fs:(1/fs)*(length(POz)-1));

% Generate sound
sineWave = sin(2*pi*50*time);
% Generate noise
v = 0.005; % variance
noise = sqrt(v)*randn(1,length(POz));
% create noisy sine wave
ref = sineWave + noise;

% define learning rates and model orders to be tested
mus = [0.001,0.01,0.1];
Ms = [5,10,20];

windowLength = 6144;
overlap = 2/3;
nfft = 5*windowLength;

figure
spectrogram(POz, hanning(windowLength), round(overlap * windowLength), nfft, fs, 'yaxis');
title('EEG: POz Noise-Corrupted Spectrogram','fontsize',15);
ax = gca;
ax.FontSize = 15;
xlabel('Time (mins)')
ylabel('Frequency (Hz)')
% yticks(0:10:60);
ylim([0 60]);
c = colorbar('TickLabelInterpreter', 'latex');
c.Label.String = "Power (dB)";
c.Label.Interpreter = 'latex';
set(gcf,'color','w');

%% iterating over learning rates and model orders, using ANC-LMS
count = 1;
figure
hold on
for i=1:length(Ms)
    for j=1:length(mus)
        % ANC-LMS
        [~,xhatANC,~] = LMS_ANC(POz',ref,mus(j),Ms(i),0); % (x_input,z_input,mu,order,leakage,delay)
         
        subplot(3,3,count) % 3 by 3 subplot
        spectrogram(xhatANC, hanning(windowLength), round(overlap * windowLength), nfft, fs, 'yaxis');
        title(sprintf('POz Spectrogram: $\\mathbf{M}$ = %0.0f and $\\mathbf{\\mu}$ = %0.3f',Ms(i),mus(j)),'Interpreter','latex','fontsize',12)
        ax = gca;
        ax.FontSize = 13;
        xlabel('Time (mins)','fontsize',13)
        ylabel('Frequency (Hz)','fontsize',13)
        yticks(0:10:60);
        ylim([0 60]);
        c = colorbar('TickLabelInterpreter', 'latex');
        c.Label.String = "Power (dB)";
        c.Label.Interpreter = 'latex';
        set(gcf,'color','w');

        disp(['Plot: ',num2str(count)]);
        count = count + 1;
    end
end


