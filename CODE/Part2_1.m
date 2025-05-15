clear all
close all
clc;

addpath Algorithms/

f = FunctionCollection();

%% Part(b)

% Initialisation=
% Number of samples
N = 1000;
% Autoregressive coefficients
% x(n) - a1x(n-1) - a2x(n-2) = noise
a0 = 1;
a1 = -0.1;
a2 = -0.8;
coeff = [a0 a1 a2];
% Variance
noiseVar = 0.25;
% Number of trials
nTrials = 100;

% % noise and x setup within loop
% wgn = sqrt(noiseVar).*randn(nTrials,N);
% data =  filter(1,coeff,wgn);

% LMS steps
step_sizes = [0.01, 0.05];
% Number of LMS steps
nSteps = length(step_sizes);

colour = {'b','r','g'};

figure;
subplot(1,2,1);
hold on
% LMS algorithm
for step = 1:nSteps
    % Generate data for current trial
    wgn = sqrt(noiseVar).*randn(1,N);
    data =  filter(1,coeff,wgn);
    
    % Apply LMS algorithm
    [y_hat,error,~] = LMS(data,data,step_sizes(step),length(coeff)-1,0);
    plot(10*log10(error.^2),'color',colour{step},'LineWidth',1.5)
end
grid minor
xlabel('Sample')
ylabel('Squared error (dB)')
legend('$\mu$ = 0.01', '$\mu$ = 0.05', 'Interpreter', 'latex')
title('Squared Prediction Error of LMS: 1 realisation')

% Error storage
avg_error = zeros(nSteps, N);
sqr_error = zeros(nTrials, N); % error.^2
misAd = zeros(nSteps,1);
weights = zeros(length(coeff)-1,N);
avg_weights = zeros(2*(length(coeff)-1),N);

% LMS algorithm
for step = 1:nSteps
    weights = zeros(length(coeff)-1,N);
    for trial = 1:nTrials
        % Generate data for current trial
        wgn = sqrt(noiseVar).*randn(1,N);
        data =  filter(1,coeff,wgn);
        
        % Apply LMS algorithm
        [y_hat,error,w] = LMS(data,data,step_sizes(step),length(coeff)-1,0);
        sqr_error(trial,:) = error;
        weights = weights + w;
    end
    % Calculate mean squared error (MSE)
    avg_error(step, :) = mean(10*log10(sqr_error.^2),1);

    % Calculate Misadjustment
    misAd(step) = mean(mean(sqr_error(:,500:end).^2,2)/0.25-1);

    % % Calculate average weights
    avg_weights(step*2-1:step*2,:) = weights/nTrials;
end

subplot(1,2,2)
for step_idx = 1:nSteps
    plot(avg_error(step_idx,:), 'LineWidth', 2, Color=colour{step_idx})
    hold on
end
grid minor
xlabel('Sample')
ylabel('MSE (dB)')
legend('$\mu$ = 0.01', '$\mu$ = 0.05', 'Interpreter', 'latex')
title('Learning Curve of LMS: 100 realisations')

disp(['Misadjustment values: timestep-',num2str(step_sizes(1)),'=',num2str(misAd(1)), ...
    '   timestep-',num2str(step_sizes(2)),'=',num2str(misAd(2))]);

%% Part(d)
figure;
subplot(1,2,1)
plot(avg_weights(1,:), 'LineWidth', 2, Color='b',DisplayName='a1')
hold on
plot(avg_weights(2,:), 'LineWidth', 2, Color='r',DisplayName='a2')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{1}$ and $\\hat{a}_{2}$}: 100 Realisations, $\\mu$ = %0.2f',0.01),'interpreter','latex','fontsize',14);
legend;

subplot(1,2,2)
plot(avg_weights(3,:), 'LineWidth', 2, Color='b',DisplayName='a1')
hold on
plot(avg_weights(4,:), 'LineWidth', 2, Color='r',DisplayName='a2')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{1}$ and $\\hat{a}_{2}$}: 100 Realisations, $\\mu$ = %0.2f',0.05),'interpreter','latex','fontsize',14);
legend;


%% Part(f) leaky mu=0.01


% Initialisation=
% Number of samples
N = 1000;
% Autoregressive coefficients
% x(n) - a1x(n-1) - a2x(n-2) = noise
a0 = 1;
a1 = -0.1;
a2 = -0.8;
coeff = [a0 a1 a2];
% Variance
noiseVar = 0.25;
% Number of trials
nTrials = 100;

% % noise and x setup within loop
% wgn = sqrt(noiseVar).*randn(nTrials,N);
% data =  filter(1,coeff,wgn);

% LMS steps
step_size = 0.01;
leakage = [0.1 0.5 0.9];
% Number of LMS steps
nSteps = length(step_sizes);

colour = {'b','r','g'};

% Error storage
avg_error = zeros(nSteps, N);
sqr_error = zeros(nTrials, N); % error.^2
misAd = zeros(nSteps,1);
weights = zeros(length(coeff)-1,N);
avg_weights = zeros(2*(length(coeff)-1),N);

% LMS algorithm
for l = 1:length(leakage)
    weights = zeros(length(coeff)-1,N);
    for trial = 1:nTrials
        % Generate data for current trial
        wgn = sqrt(noiseVar).*randn(1,N);
        data =  filter(1,coeff,wgn);
        
        % Apply LMS algorithm
        [y_hat,error,w] = LMS(data,data,step_size,length(coeff)-1,leakage(l));
        weights = weights + w;
    end

    % % Calculate average weights
    avg_weights(l*2-1:l*2,:) = weights/nTrials;
end
 
figure;
subplot(1,2,1)
plot(avg_weights(1,:), 'LineWidth', 2, Color='b',DisplayName='\mu=0.01, \gamma=0.1')
hold on
plot(avg_weights(3,:), 'LineWidth', 2, Color='r',DisplayName='\mu=0.01, \gamma=0.5')
hold on
plot(avg_weights(5,:), 'LineWidth', 2, Color='g',DisplayName='\mu=0.01, \gamma=0.9')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{1}$}: 100 Realisations'),'interpreter','latex','fontsize',14);
legend(Location="southeast");

subplot(1,2,2)
plot(avg_weights(2,:), 'LineWidth', 2, Color='b',DisplayName='\mu=0.01, \gamma=0.1')
hold on
plot(avg_weights(4,:), 'LineWidth', 2, Color='r',DisplayName='\mu=0.01, \gamma=0.5')
hold on
plot(avg_weights(6,:), 'LineWidth', 2, Color='g',DisplayName='\mu=0.01, \gamma=0.9')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{2}$}: 100 Realisations'),'interpreter','latex','fontsize',14);
legend(Location="southeast");

%% Part(f) leaky mu=0.05

% Initialisation=
% Number of samples
N = 1000;
% Autoregressive coefficients
% x(n) - a1x(n-1) - a2x(n-2) = noise
a0 = 1;
a1 = -0.1;
a2 = -0.8;
coeff = [a0 a1 a2];
% Variance
noiseVar = 0.25;
% Number of trials
nTrials = 100;

% % noise and x setup within loop
% wgn = sqrt(noiseVar).*randn(nTrials,N);
% data =  filter(1,coeff,wgn);

% LMS steps
step_size = 0.05;
leakage = [0.1 0.5 0.9];
% Number of LMS steps
nSteps = length(step_sizes);

colour = {'b','r','g'};

% Error storage
avg_error = zeros(nSteps, N);
sqr_error = zeros(nTrials, N); % error.^2
misAd = zeros(nSteps,1);
weights = zeros(length(coeff)-1,N);
avg_weights = zeros(2*(length(coeff)-1),N);

% LMS algorithm
for l = 1:length(leakage)
    weights = zeros(length(coeff)-1,N);
    for trial = 1:nTrials
        % Generate data for current trial
        wgn = sqrt(noiseVar).*randn(1,N);
        data =  filter(1,coeff,wgn);
        
        % Apply LMS algorithm
        [y_hat,error,w] = LMS(data,data,step_size,length(coeff)-1,leakage(l));
        weights = weights + w;
    end

    % % Calculate average weights
    avg_weights(l*2-1:l*2,:) = weights/nTrials;
end
 
figure;
subplot(1,2,1)
plot(avg_weights(1,:), 'LineWidth', 2, Color='b',DisplayName='\mu=0.05, \gamma=0.1')
hold on
plot(avg_weights(3,:), 'LineWidth', 2, Color='r',DisplayName='\mu=0.05, \gamma=0.5')
hold on
plot(avg_weights(5,:), 'LineWidth', 2, Color='g',DisplayName='\mu=0.05, \gamma=0.9')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{1}$}: 100 Realisations'),'interpreter','latex','fontsize',14);
legend(Location="southeast");

subplot(1,2,2)
plot(avg_weights(2,:), 'LineWidth', 2, Color='b',DisplayName='\mu=0.05, \gamma=0.1')
hold on
plot(avg_weights(4,:), 'LineWidth', 2, Color='r',DisplayName='\mu=0.05, \gamma=0.5')
hold on
plot(avg_weights(6,:), 'LineWidth', 2, Color='g',DisplayName='\mu=0.05, \gamma=0.9')
grid minor
xlabel('Sample')
ylabel('Weight values')
title(sprintf('\\textbf{Time Evolution of $\\hat{a}_{2}$}: 100 Realisations'),'interpreter','latex','fontsize',14);
legend(Location="southeast");
