clear all
close all
clc;

addpath Algorithms/

%%
load('Data/time-series.mat')

N = length(y);
mu = 0.00001;
gamma = 0;
order = 4;
y = y -mean(y);

[yhat,error,~] = LMS(y',y',mu,order,gamma);

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));
disp(['MSE= ',num2str(MSE_db),'dB    R_p = ',num2str(R_p)]);

MSE_end_db = mean(error(800:1000).^2);
MSE_end_db = 10*log10(MSE_end_db);
R_end_p = 10*log10(var(yhat(800:1000))/var(error(800:1000)));
disp(['For the end portion: MSE= ',num2str(MSE_end_db),'dB    R_p = ',num2str(R_end_p)]);

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));
disp(['MSE= ',num2str(MSE_db),'dB    R_p = ',num2str(R_p)]);

MSE_end_db = mean(error(800:1000).^2);
MSE_end_db = 10*log10(MSE_end_db);
R_end_p = 10*log10(var(yhat(800:1000))/var(error(800:1000)));
disp(['For the end portion: MSE= ',num2str(MSE_end_db),'dB    R_p = ',num2str(R_end_p)]);


figure
subplot(1,2,1)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('MSE=%.2f dB, Rp=%.2f dB',MSE_db,R_p),'fontsize',8,'interpreter','latex')
legend
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('End: MSE=%.2f dB, Rp=%.2f dB',MSE_end_db,R_end_p),'fontsize',8,'interpreter','latex')
legend
sgtitle('One-Step-Ahead Prediction of AR(4) Time Series','fontsize',20)
xlim([800 1000])
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([800 1000])


%% Dynamical Perceptron
clear
load('Data/time-series.mat')

mu = 0.00001;
order = 4;
y = y -mean(y);
w_init = zeros(order,1);
bias = false;
alpha = 1;

[yhat,error,~,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha);

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));
disp(['MSE= ',num2str(MSE_db),'dB    R_p = ',num2str(R_p)]);

MSE_end_db = mean(error(800:1000).^2);
MSE_end_db = 10*log10(MSE_end_db);
R_end_p = 10*log10(var(yhat(800:1000))/var(error(800:1000)));
disp(['For the end portion: MSE= ',num2str(MSE_end_db),'dB    R_p = ',num2str(R_end_p)]);


figure
subplot(1,2,1)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('MSE=%.2f dB, Rp=%.2f dB',MSE_db,R_p),'fontsize',8,'interpreter','latex')
legend
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('End: MSE=%.2f dB, Rp=%.2f dB',MSE_end_db,R_end_p),'fontsize',8,'interpreter','latex')
legend
sgtitle('One-Step-Ahead Prediction of AR(4) Time Series','fontsize',20)
xlim([800 1000])
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([800 1000])

%% Scaling the alpha
clear
load('Data/time-series.mat')

mu = 0.0000001;
order = 4;
y = y -mean(y);
w_init = zeros(order,1);
bias = false;
alpha = 50:0.1:100;

MSEs = zeros(1,length(alpha));
R_ps = MSEs;

for i = 1:length(alpha)
    [yhat,error,~,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha(i));

    % calculating the MSE between the true and estimated signals
    MSEs(i) = 10*log10(mean(abs(error).^2));
    R_ps(i) = 10*log10(var(yhat)/var(error));
end


figure
subplot(1,2,1)
plot(alpha,MSEs,'b','LineWidth',1.5)
hold on
[val,ind] = min(MSEs);
plot(alpha(1) + ind*0.1,val,'r*','MarkerSize',20)
xlabel('alpha')
ylabel('Mean Squared Error (dB)')
title(sprintf('Min: $\\alpha$= %.1f',alpha(1) + ind*0.1),'fontsize',15,'interpreter','latex')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([alpha(1) alpha(end)])

subplot(1,2,2)
plot(alpha,R_ps,'b','LineWidth',1.5)
hold on
[val,ind] = max(R_ps);
plot(alpha(1) + ind*0.1,val,'r*','MarkerSize',20)
xlabel('alpha')
ylabel('Prediction Gain (dB)')
title(sprintf('Max: $\\alpha$= %.1f',alpha(1) + ind*0.1),'fontsize',15,'interpreter','latex')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([alpha(1) alpha(end)])

sgtitle(sprintf('Finding the Optimal Value for $\\alpha$, $\\mu$ = %.1e',mu),'fontsize',20,'interpreter','latex')

%% plot with optimal parameters

mu = 0.0000001;
order = 4;
y = y -mean(y);
w_init = zeros(order,1);
bias = false;
alpha = 92;

[yhat,error,w,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha);

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));
disp(['MSE= ',num2str(MSE_db),'dB    R_p = ',num2str(R_p)]);

MSE_end_db = mean(error(750:1000).^2);
MSE_end_db = 10*log10(MSE_end_db);
R_end_p = 10*log10(var(yhat(750:1000))/var(error(750:1000)));
disp(['For the end portion: MSE= ',num2str(MSE_end_db),'dB    R_p = ',num2str(R_end_p)]);

figure
subplot(1,2,1)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
title(sprintf('MSE=%.2f dB, Rp=%.2f dB',MSE_db,R_p),'fontsize',8,'interpreter','latex')
xlabel('Time (n)')
ylabel('(AU)')
legend
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
title(sprintf('End: MSE=%.2f dB, Rp=%.2f dB',MSE_end_db,R_end_p),'fontsize',8,'interpreter','latex')
xlabel('Time (n)')
ylabel('(AU)')
legend
xlim([800 1000])
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([800 1000])

sgtitle('One-Step-Ahead Prediction of AR(4) Time Series','fontsize',20)

% plot the evolution of the weights
figure
plot(w(1,:),'b','LineWidth',1.5,DisplayName='w1')
hold on
plot(w(2,:),'r','LineWidth',1.5,DisplayName='w2')
hold on
plot(w(3,:),'m','LineWidth',1.5,DisplayName='w3')
hold on
plot(w(4,:),'c','LineWidth',1.5,DisplayName='w4')
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('Evolution of LMS-DP Weights (Scaled)','Interpreter','latex','fontsize',18)
legend(Location='northwest')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor

%% Add bias
clear
load('Data/time-series.mat')

mu = 0.0000001;
order = 4;
y = y -mean(y);
bias = true;
biasDC= 10;
w_init = zeros(order+bias,1);
alpha = 92;

[yhat,error,w,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha,biasDC);

% calculating the MSE between the true and estimated signals
MSE = mean(error.^2);
MSE_db = 10*log10(MSE);
R_p = 10*log10(var(yhat)/var(error));
disp(['MSE= ',num2str(MSE_db),'dB    R_p = ',num2str(R_p)]);

MSE_end_db = mean(error(750:1000).^2);
MSE_end_db = 10*log10(MSE_end_db);
R_end_p = 10*log10(var(yhat(750:1000))/var(error(750:1000)));
disp(['For the end portion: MSE= ',num2str(MSE_end_db),'dB    R_p = ',num2str(R_end_p)]);

figure
subplot(1,2,1)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
title(sprintf('MSE=%.2f dB, Rp=%.2f dB',MSE_db,R_p),'fontsize',8,'interpreter','latex')
xlabel('Time (n)')
ylabel('(AU)')
legend
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
title(sprintf('End: MSE=%.2f dB, Rp=%.2f dB',MSE_end_db,R_end_p),'fontsize',8,'interpreter','latex')
xlabel('Time (n)')
ylabel('(AU)')
legend
xlim([800 1000])
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([800 1000])

sgtitle(sprintf('One-Step-Ahead Prediction of AR(4) Time Series: bias=%d',biasDC),'fontsize',20)

% plot the evolution of the weights
figure
plot(w(1,:),'b','LineWidth',1.5,DisplayName='w0')
hold on
plot(w(2,:),'r','LineWidth',1.5,DisplayName='w1')
hold on
plot(w(3,:),'m','LineWidth',1.5,DisplayName='w2')
hold on
plot(w(4,:),'c','LineWidth',1.5,DisplayName='w3')
hold on
plot(w(5,:),'g','LineWidth',1.5,DisplayName='w4')
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('Evolution of LMS-DP Weights (Scaled, Biased)','Interpreter','latex','fontsize',18)
legend(Location='northwest')
ax = gca;
ax.FontSize = 15; 
grid on
grid minor

%% pre-train the weights

clear
load('Data/time-series.mat')

mu = 0.0000001;
order = 4;
y = y -mean(y);
bias = true;
w_init = zeros(order+bias,1);
alpha_list = 50:0.1:90;

epochs = 100;
segLen = 20;
all_winit = zeros(order+bias,length(alpha_list));


MSEs = zeros(1,length(alpha_list));
R_ps = MSEs;

for a = 1:length(alpha_list)
    disp(['alpha=',num2str(alpha_list(a))])
    for e = 1:epochs
        y_samp = y(1:segLen);  % get sample

        % pretraining the weights
        y_hat = zeros(length(y_samp),1);
        error = y_hat;

        w = zeros(order+bias,length(y_samp)+1);
        w(:,1) = w_init;

        x_delayed = zeros(order,length(y_samp));
        % creating two shifted vectors by order length
        for i = 1: order
            x_delayed(i,:) = [ zeros(1,i), y_samp(1: length(y_samp)-i)']; 
        end
        if bias
            x_delayed = [ones(1,length(y_samp)); x_delayed];
        end

        for k = 1: length(y_samp)
            % calculate the prediction
            y_hat(k) = alpha_list(a)*tanh(w(:,k)'*x_delayed(:,k));
            error(k) = y_samp(k)-y_hat(k);
            act_function = alpha_list(a)*(1-(y_hat(k)./alpha_list(a))^2);

            % update
            w(:,k+1)=w(:,k)+(mu*act_function*error(k)).*x_delayed(:,k);
        end
        w_init = w(:,end);
    end
    all_winit(:,a) = w_init;
    [yhat,err,~,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha_list(a));

    % calculating the MSE between the true and estimated signals
    MSEs(a) = 10*log10(mean(abs(err(order+1:end)).^2));
    R_ps(a) = 10*log10(var(yhat(order+1:end))/var(err(order+1:end)));

end

figure
subplot(1,2,1)
plot(alpha_list,MSEs,'b','LineWidth',1.5)
hold on
[val,ind] = min(MSEs);
plot(alpha_list(1) + ind*0.1,val,'r*','MarkerSize',20)
xlabel('alpha')
ylabel('Mean Squared Error (dB)')
title(sprintf('Min: $\\alpha$= %.1f',alpha_list(1) + ind*0.1),'fontsize',15,'interpreter','latex')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([alpha_list(1) alpha_list(end)])

subplot(1,2,2)
plot(alpha_list,R_ps,'b','LineWidth',1.5)
hold on
[val,ind] = max(R_ps);
plot(alpha_list(1) + ind*0.1,val,'r*','MarkerSize',20)
xlabel('alpha')
ylabel('Prediction Gain (dB)')
title(sprintf('Max: $\\alpha$= %.1f',alpha_list(1) + ind*0.1),'fontsize',15,'interpreter','latex')
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([alpha_list(1) alpha_list(end)])

sgtitle(sprintf('Finding the Optimal Value for $\\alpha$, $\\mu$ = %.1e',mu),'fontsize',20,'interpreter','latex')
%%

N = length(y);
n = 1: N;
mu = 0.0000001;
gamma = 0;
order = 4;
y = y -mean(y);
bias = true;
alpha = 70;
w_init = all_winit(:,alpha_list == alpha);

[y_hat,error,w,~] = LMS_DP(y,y,mu,order,w_init,bias,alpha);

% calculating the MSE between the true and estimated signals
MSE = 10*log10(mean(abs(error).^2));
MSE_db = 10*log10(MSE);

MSE_end_db = mean(mean(abs(error(800:1000)).^2));
MSE_end_db = 10*log10(MSE_end_db);

R_p = 10*log10(var(y_hat)/var(error));
R_end_p = 10*log10(var(y_hat(800:1000))/var(error(800:1000)));

figure
subplot(1,2,1)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('MSE=%.2f dB, Rp=%.2f dB',MSE_db,R_p),'fontsize',8,'interpreter','latex')
legend
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')

subplot(1,2,2)
plot(y,'b','LineWidth',1.5,DisplayName='Centred Time Series')
hold on
plot(yhat,'r','LineWidth',1.5,DisplayName='Estimated Time Series')
xlabel('Time (n)')
ylabel('(AU)')
title(sprintf('End: MSE=%.2f dB, Rp=%.2f dB',MSE_end_db,R_end_p),'fontsize',8,'interpreter','latex')
legend
sgtitle('One-Step-Ahead Prediction of AR(4) Time Series','fontsize',20)
xlim([800 1000])
ax = gca;
ax.FontSize = 18; 
grid on
grid minor
set(gcf,'color','w')
xlim([800 1000])

% plot the evolution of the weights
figure
plot(w(1,:),'b','LineWidth',1.5,DisplayName='w0')
hold on
plot(w(2,:),'r','LineWidth',1.5,DisplayName='w1')
hold on
plot(w(3,:),'m','LineWidth',1.5,DisplayName='w2')
hold on
plot(w(4,:),'c','LineWidth',1.5,DisplayName='w3')
hold on
plot(w(5,:),'g','LineWidth',1.5,DisplayName='w4')
xlabel('Time Index (n)','fontsize',18)
ylabel('(AU)','fontsize',18)
title('Evolution of LMS-DP Weights (Scaled, Biased)','Interpreter','latex','fontsize',18)
legend(Location='northwest')
ylim([-0.06 0.06])
ax = gca;
ax.FontSize = 15; 
grid on
grid minor