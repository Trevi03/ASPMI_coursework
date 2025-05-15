clear all
close all
clc;

%% Part (a)
% [xRRI,fsRRI]=ECG_to_RRI(data1,500);
% data = readtable("../Group_1_2025-01-21-160844_EEG.xlsx");
% extract_data = data.CH1;

load Data/RRI_data.mat


xRRI1 = detrend(xRRI_trial1);
xRRI2 = detrend(xRRI_trial2);
xRRI3 = detrend(xRRI_trial3);

xRRI = {xRRI1;xRRI2;xRRI3};
fsRRI = 4;

for i = 1:3
    xRRIi = xRRI{i};
    len = length(xRRIi);
    [psdStan,fStan] = pwelch(xRRIi, hamming(len), 0,1024, 4, 'onesided');

    [psd150,f150] = pwelch(xRRIi, hamming(150*fsRRI), 0,1024, 4, 'onesided');
    [psd50,f50] = pwelch(xRRIi, rectwin(50*fsRRI), 0,1024, 4, 'onesided');
    
    subplot(3,1,i)  % on the right column
    plot(fStan,10*log10(psdStan),'b','LineWidth',2) 
    hold on
    plot(f150,10*log10(psd150),'g','LineWidth',2)
    hold on
    plot(f50,10*log10(psd50),'r','LineWidth',2)
    ax = gca;
    ax.FontSize = 12;
    ylim([-100 0]);
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB)')
    title(sprintf('Averaged Periodogram: RRI%d',i),'fontsize',15)
    grid on
    grid minor
end
legend('Standard','150 s window','50 s window')
set(gcf,'color','w')


%% Part (c)

order = [4,10,7];
RRI = {xRRI1,xRRI2,xRRI3}; % [4, ,6]

figure; hold on
for i=1:length(order)
    [stanPSD,stanf] = pwelch(RRI{i}, length(RRI{i}), 0,1024, 4, 'onesided');
    [psdEst, w] = pyulear(RRI{i}, order(i), 2048, 4);
    subplot(1,3,i)
    plot(stanf, 10*log10(stanPSD),LineWidth=2);
    hold on
    plot(w, 10*log10(psdEst),LineWidth=3,Color='r');
    grid on
    grid minor
    set(gca,'fontsize', 14);
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB)'); 
    title(['RRI',num2str(i),' Data: AR(',num2str(order(i)),')']);
end

legend('Standard', 'AR(i)'); hold off;
% legend;


