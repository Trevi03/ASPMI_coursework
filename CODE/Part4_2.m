clear all
close all
clc;

addpath Algorithms/
%%

t = ((2*pi)/100):((2*pi)/100):10*pi;
y1 = sin(t);
y2 = sin(0.5*t);
y3 = sin(4*t);

y = y1 + rand(1,length(t));

% Compute convolutions
auto_conv_y1 = conv(y, y1, 'same'); % Auto-convolution of y1
conv_y1_y2 = conv(y, y2, 'same');   % Convolution of y1 with y2
conv_y1_y3 = conv(y, y3, 'same');   % Convolution of y1 with y3

time = 1:length(y1);

% Plot the results
figure;
plot(time,auto_conv_y1,LineWidth=2.0,Color="#D95319",DisplayName='auto y1');
hold on
plot(time,conv_y1_y2,LineWidth=2.0,Color="#0072BD",DisplayName='y1 with y2');
hold on
plot(time,conv_y1_y3,LineWidth=2.0,Color="#77AC30",DisplayName='y1 with y3');
title('Convolution Plots');
xlabel('Time');
ylabel('Amplitude');
legend
grid minor

%% Signal identification

% change as needed
[signal_class,conv_y_y1,conv_y_y2,conv_y_y3,time] = sigID(y1);

figure;
subplot(3,1,1);
plot(time,conv_y_y1,LineWidth=2.0,Color="#D95319",DisplayName='y and y1');
hold on
plot(time,conv_y_y2,LineWidth=2.0,Color="#0072BD",DisplayName='y and y2');
hold on
plot(time,conv_y_y3,LineWidth=2.0,Color="#77AC30",DisplayName='y and y3');
title('y=y1 system identification');
title(sprintf('y=y1 system identification: detected signal = y%d',signal_class),'fontsize',15,'interpreter','latex')
xlabel('Time');
ylabel('Amplitude');
legend
xlim([time(1) time(end)])
grid minor

% change as needed
[signal_class,conv_y_y1,conv_y_y2,conv_y_y3,time] = sigID(y2);

subplot(3,1,2);
plot(time,conv_y_y1,LineWidth=2.0,Color="#D95319",DisplayName='y and y1');
hold on
plot(time,conv_y_y2,LineWidth=2.0,Color="#0072BD",DisplayName='y and y2');
hold on
plot(time,conv_y_y3,LineWidth=2.0,Color="#77AC30",DisplayName='y and y3');
title('y=y1 system identification');
title(sprintf('y=y2 system identification: detected signal = y%d',signal_class),'fontsize',15,'interpreter','latex')
xlabel('Time');
ylabel('Amplitude');
legend
xlim([time(1) time(end)])
grid minor

% change as needed
[signal_class,conv_y_y1,conv_y_y2,conv_y_y3,time] = sigID(y3);

subplot(3,1,3);
plot(time,conv_y_y1,LineWidth=2.0,Color="#D95319",DisplayName='y and y1');
hold on
plot(time,conv_y_y2,LineWidth=2.0,Color="#0072BD",DisplayName='y and y2');
hold on
plot(time,conv_y_y3,LineWidth=2.0,Color="#77AC30",DisplayName='y and y3');
title('y=y1 system identification');
title(sprintf('y=y3 system identification: detected signal = y%d',signal_class),'fontsize',15,'interpreter','latex')
xlabel('Time');
ylabel('Amplitude');
legend
xlim([time(1) time(end)])
grid minor



function [signal_class,conv_y_y1,conv_y_y2,conv_y_y3,t] = sigID(y)
    % initialise signals
    t = ((2*pi)/100):((2*pi)/100):10*pi;
    y1 = sin(t);
    y2 = sin(0.5*t);
    y3 = sin(4*t);

    % Flip reference signals for cross-correlation (flipped convolution)
    y1_flipped = flip(y1);
    y2_flipped = flip(y2);
    y3_flipped = flip(y3);
    
    % Perform convolution with flipped signals
    conv_y_y1 = conv(y, y1_flipped, 'same');
    conv_y_y2 = conv(y, y2_flipped, 'same');
    conv_y_y3 = conv(y, y3_flipped, 'same');
    
    % Apply absolute value to remove negative correlations
    corr_y_y1 = abs(conv_y_y1);
    corr_y_y2 = abs(conv_y_y2);
    corr_y_y3 = abs(conv_y_y3);
    
    % Compute maximum correlation values
    max_corr_y1 = max(corr_y_y1);
    max_corr_y2 = max(corr_y_y2);
    max_corr_y3 = max(corr_y_y3);
    
    % Identify which signal y is closest to
    [~, signal_class] = max([max_corr_y1, max_corr_y2, max_corr_y3]);
end


