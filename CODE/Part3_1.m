clear all
close all
clc;

addpath Algorithms/

%% Part(a)

N = 1000;
nTrials = 100;
mu = 0.01; 
b0 = 1;
b = [complex(1.5,1), complex(2.5,-0.5)];
order = length(b); % Order

% CLM coefficinets, error and signal
h = complex(zeros(nTrials, order, N)); 
e = complex(zeros(nTrials, N));
y = complex(zeros(1,N));

% ACLM coefficinets, error and signal
ah = h;
ag = h;
ae = e;

for k=1:nTrials
    % WLMA(1)
    x = randn(1,N)+1j*randn(1,N); 
    y(2:N) = b0*x(2:N) + b(1)*x(1:N-1) + b(2)*conj(x(1:N-1));
    % y(n) = x(n) + b1x(n−1) + b2x^*(n−1) x∼N(0,1)
    
    % implement clms algorithm 
    [h(k,:,:),e(k,:)] = CLMS(x,y,mu,order);
    % implement aclms algorithm
    [ah(k,:,:),ag(k,:,:),ae(k,:)] = ACLMS(x,y,mu,order);
end

figure
hold on
subplot(1,2,1)
scatter(real(y),imag(y),20,'filled')
hold on
scatter(real(x),imag(x),20,'filled','r')
xlabel("Re",'fontsize',18)
ylabel("Im",'fontsize',18)
title('Circularity','fontsize',20)
xlim([-10 10])
ylim([-10 10])
legend('WLMA(1)','WGN')
ax = gca;
ax.FontSize = 15;
grid on
grid minor
set(gcf,'color','w')

error = abs(e).^2;
error = squeeze(mean(error, 1));
errorA = abs(ae).^2;
errorA = squeeze(mean(errorA, 1));

subplot(1,2,2)
plot(10*log10(error),'LineWidth',2)
hold on
plot(10*log10(errorA),'r','LineWidth',2)
title('Learning Curves CLMS and ACLMS Algorithms', 'fontsize', 20)
xlabel('Sample (n)')
ylabel('Squared Error(dB)')
legend('CLMS','ACLMS')
ax = gca;
ax.FontSize = 15;
grid on
grid minor

%% Part (b)
clear

lowWind = load('wind-dataset/high-wind.mat');
medWind = load('wind-dataset/medium-wind.mat');
highWind = load('wind-dataset/low-wind.mat');

% represent north and south v as complex array
v_low = complex(lowWind.v_east,lowWind.v_north);
v_med = complex(medWind.v_east,medWind.v_north);
v_high = complex(highWind.v_east,highWind.v_north);

% find out the circularity
circ_low = abs(mean((v_low).^2)/mean(abs(v_low).^2));
circ_med = abs(mean((v_med).^2)/mean(abs(v_med).^2));
circ_high =  abs(mean((v_high).^2)/mean(abs(v_high).^2));

disp('Circularity:')
disp(['Low speed = ',num2str(circ_low)])
disp(['Medium speed = ',num2str(circ_med)])
disp(['High speed = ',num2str(circ_high)])

% centre of mass
com_low = [mean(real(v_low)),mean(imag(v_low))];
com_med = [mean(real(v_med)),mean(imag(v_med))];
com_high = [mean(real(v_high)),mean(imag(v_high))];

figure;
% low
subplot(1,3,1)
scatter(real(v_low),imag(v_low),10,'filled','r')
hold on
scatter(com_low(1),com_low(2),30,'filled','y')
xlabel("Re (East)")
ylabel("Im (North)")
title(sprintf('Low-speed Data: $\\rho$ = %0.2f',circ_low),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor

% med
subplot(1,3,2)
scatter(real(v_med),imag(v_med),10,'filled','m')
hold on
scatter(com_med(1),com_med(2),30,'filled','y')
xlabel("Re (East)")
ylabel("Im (North)")
title(sprintf('Medium-speed Data: $\\rho$ = %0.2f',circ_med),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor

% high
subplot(1,3,3)
scatter(real(v_high),imag(v_high),10,'filled','b')
hold on
scatter(com_high(1),com_high(2),30,'filled','y')
xlabel("Re (East)")
ylabel("Im (North)")
title(sprintf('High-speed Data: $\\rho$ = %0.2f',circ_high),'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16;
grid on
grid minor

%% Configure the CLMS and ACLMS filters in a prediction setting

N = length(v_low);
nTrials = 100;

mu = [0.001 0.005 0.1]; 
filtLen = 25;

all_v = [v_low, v_med, v_high];

for k=1:3 % for 3 speeds
    e = complex(zeros(length(filtLen), N));
    ae = e;
    y = all_v(:,k);

    for f=1:filtLen
        x = [0; y(1:end-1)]'; % time shifted x in prediction mode for predicting future values
        
        % implement clms algorithm 
        [~,e(f,:)] = predCLMS(x,y',mu(k),f);
        % implement aclms algorithm
        [~,~,ae(f,:)] = predACLMS(x,y',mu(k),f);
    end

    error = abs(e).^2;
    error = nanmean(error, 2);
    errorA = abs(ae).^2;
    errorA = nanmean(errorA, 2);

    subplot(1,3,k)
    plot(10*log10(error),'LineWidth',2,DisplayName="CLMS")
    [~,idx1] = min(10*log10(error));
    hold on
    plot(10*log10(errorA),'r','LineWidth',2,DisplayName="ACLMS")
    [~,idx2] = min(10*log10(errorA));    
    if k == 1
    	title('Learning Curves - Low Speed', 'fontsize', 16)
        disp('Learning Curves - Low Speed')
    elseif k == 2
        title('Learning Curves - Mid Speed', 'fontsize', 16)
        disp('Learning Curves - Mid Speed')
    else
        title('Learning Curves - High Speed', 'fontsize', 16)
        disp('Learning Curves - High Speed')
    end
    disp(['Minimum indexes: CLMS=',num2str(idx1),'  ACLMS=',num2str(idx2)]);

    xlabel('Filter Order')
    ylabel('MSPE (dB)')
    ax = gca;
    ax.FontSize = 15;
    legend;
    grid on
    grid minor
    axis tight
end


