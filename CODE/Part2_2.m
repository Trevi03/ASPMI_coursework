clear all
close all
clc;

addpath Algorithms/

%% Part (a)
N =  1000;
nTrials = 1000;
rho = 0.001;
mu = [0.01,0.1,0.2,0.2,0.2];
alpha = 0.8;

% define the MA model
b = [1,0.9];
order = length(b) -1;
a = 1;
v = 0.5;

% initiate noise
eta = zeros(nTrials,N);
for i = 1:nTrials
    eta(i,:) = sqrt(v)*randn(1,N);
end

figure;
hold on
colors = {'b','r','c','m',"#7E2F8E"};
subplot(1,3,1)
for i = 1:5
    w_all = zeros(nTrials,N);
    
    for n = 1: nTrials
       x = filter(b,a,eta(n,:));
       [~,~,w,~] = LMS_GASS(eta(n,:),x,mu(i),rho,alpha,order,i);
       w_all(n,:) = w; 
    end
    w_error = b(2)*ones(nTrials,N) - w_all;
    plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
    hold on
end
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',12);
xlabel('Sample (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0.2$','interpreter','latex')
ylim([-0.05 1]);
grid on
grid minor
%%
mu = [0.01,0.1,0.1,0.1,0.1];
subplot(1,3,2)
for i = 1:5
    w_all = zeros(nTrials,N);
    for n = 1: nTrials
       
       x = filter(b,a,eta(n,:));
       [xHat,error,w,musFinal] = LMS_GASS(eta(n,:),x,mu(i),rho,alpha,order,i);
       w_all(n,:) = w; 
    end
    w_error = b(2)*ones(nTrials,N) - w_all;
    plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
    hold on
end
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',13);
xlabel('Sample (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0.1$','interpreter','latex')
ylim([-0.05 1]);
grid on
grid minor
%%
mu = [0.01,0.1,0,0,0];
subplot(1,3,3)
for i = 1:5
    w_all = zeros(nTrials,N);
    for n = 1: nTrials
       
       x = filter(b,a,eta(n,:));
       [~,~,w,~] = LMS_GASS(eta(n,:),x,mu(i),rho,alpha,order,i);
       w_all(n,:) = w; 
    end
      w_error = b(2)*ones(nTrials,N) - w_all;
      plot(nanmean(w_error),'color',colors{i},'LineWidth',1.5);
      hold on
end
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',12);
xlabel('Sample (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
title('$\mu_{initGASS} = 0$','interpreter','latex')
ylim([-0.05 1]);
grid on
grid minor
sgtitle('Weight Error Curves','fontsize',16)
set(gcf,'color','w')

%%
mu = [0.01,0.1,0.1,0.1,0.1];
figure;
for i = 1:5
    w_all = zeros(nTrials,N);
    for n = 1: nTrials
       
       x = filter(b,a,eta(n,:));
       [xHat,error,w,musFinal] = LMS_GASS(eta(n,:),x,mu(i),rho,alpha,order,i);
       w_all(n,:) = w; 
    end
      w_error = b(2)*ones(nTrials,N) - w_all;
      plot(10*log10(mean(w_error.^2)),'color',colors{i},'LineWidth',1.5);
      hold on
end
legend('$\mu = 0.01$','$\mu = 0.1$','Benveniste','Ang and Farhang','Matthews and Xie','Interpreter','latex','fontsize',13);
xlabel('Sample (n)','fontsize',15)
ylabel('Squared Weight Error (dB)','fontsize',15)
title('$\mu_{initGASS} = 0.1$','interpreter','latex')
grid on
grid minor

%% Part (c)
clear all

N =  1000;
nTrials = 1000;
rho = 0.001;
alpha = 0.8;
mu = 0.1;
epsilon = 1/mu;

% define the MA model
b = [1,0.9];
order = length(b) -1;
a = 1;
v = 0.5;

% initiate noise
eta = zeros(nTrials,N);
for i = 1:nTrials
    eta(i,:) = sqrt(v)*randn(1,N);
end

N = N/2;
% initialise error and weights array for Benveniste
w_allB = zeros(nTrials,N);
w_errB = w_allB;
% initialise error and weights array for GNGD
w_allG = w_allB;
w_errG = w_allB;

for i = 1:nTrials
    % generate signal
    x = filter(b,a,eta(i,:));
    x =  x(501:end);
    
    [~,errB,wB,~] = LMS_GASS(eta(i,501:end),x,0.1,0.002,alpha,order,3);
    [~,errG,wG,~] = GNGD(eta(i,501:end),x,1,0.005,epsilon,order); 

    w_errB(i,:) = errB;
    w_allB(i,:) = wB; 
    w_allG(i,:) = wG;
    w_errG(i,:) = errG;
end

w_errorB = b(2)*ones(nTrials,N) - w_allB;
w_errorG = b(2)*ones(nTrials,N) - w_allG;


figure
subplot(1,2,1)
hold on
plot(nanmean(w_errorB),'b','LineWidth',1.5,DisplayName='Benveniste');
hold on
plot(nanmean(w_errorG),'r','LineWidth',1.5,DisplayName='GNGD');
legend('Interpreter','latex','fontsize',13);
xlabel('Sample (n)','fontsize',15)
ylabel('Weight Error (AU)','fontsize',15)
xlim([0 100])
ax = gca; ax.FontSize = 15;
title('Weight Error Curves','interpreter','latex')
grid on
grid minor

subplot(1,2,2)
plot(nanmean(10*log10(w_errorB.^2)),'b','LineWidth',1.5,DisplayName='Benveniste');
hold on
plot(nanmean(10*log10(w_errorG.^2)),'r','LineWidth',1.5,DisplayName='GNGD');
legend('Interpreter','latex','fontsize',13);
xlabel('Sample (n)','fontsize',15)
ylabel('Squared Weight Error (dB)','fontsize',15)
xlim([0 100])
ax = gca; ax.FontSize = 15;
title('Squared Weight Error Curves','interpreter','latex')
grid on
grid minor
