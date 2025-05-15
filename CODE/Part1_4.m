clear all
close all
clc;

%% Part (b)
N = 1000;
% x(n) = 2.76x(n−1)−3.81x(n−2) + 2.65x(n−3)−0.92x(n−4) + w(n)
a = [1 -2.76 3.81 -2.65 0.92];
x =  filter(1,a, randn(N,1));
x = x(500:end);

% x(isnan(x) | isinf(x)) = 0;

order = [2 4 10];  % the order with the best results
[H,w] = freqz(1,a,length(x)); % freq response
truePSD = abs(H).^2;

figure;
hold on
for i = 1: length(order)
    subplot(1,3,i);

    [arcoefs,E] = aryule(x,order(i));
    [H,w] = freqz(E^(1/2),arcoefs,length(x));
    estPSD = abs(H).^2;
    plot(w/pi,10*log10(truePSD),'b','LineWidth',2);
    hold on
    plot(w/pi,10*log10(estPSD),'r','LineWidth',2);

    tag = join(['AR(',num2str(order(i)),')']);
    legend('True',tag)
    ax.FontSize = 14;
    axis tight
    xlabel('Normalised Frequency (\pi rads/sample)')
    ylabel('PSD (dB)')
    title('AR PSD Estimation, N = 500','fontsize',15)
    grid on 
    grid minor
    set(gcf, 'color','w');
end


%% Part (c)
N = 1e4;
% x(n) = 2.76x(n−1)−3.81x(n−2) + 2.65x(n−3)−0.92x(n−4) + w(n)
a = [1 -2.76 3.81 -2.65 0.92];
x =  filter(1,a, randn(N,1));
x = x(500:end);


order = [2 4 6];  % the order with the best results
[H,w] = freqz(1,a,length(x)); % freq response
truePSD = abs(H).^2;

figure;
hold on
for i = 1: length(order)
    subplot(1,3,i);

    [arcoefs,E] = aryule(x,order(i));
    [H,w] = freqz(E^(1/2),arcoefs,length(x));
    estPSD = abs(H).^2;
    plot(w/pi,10*log10(truePSD),'b','LineWidth',2);
    hold on
    plot(w/pi,10*log10(estPSD),'r','LineWidth',2);

    tag = join(['AR(',num2str(order(i)),')']);
    legend('True',tag)
    set(gca,'fontsize', 14);
    axis tight
    xlabel('Normalised Frequency (\pi rads/sample)')
    ylabel('PSD (dB)')
    title('AR PSD Estimation, N = 9500','fontsize',15)
    grid on 
    grid minor
    set(gcf, 'color','w');
end

