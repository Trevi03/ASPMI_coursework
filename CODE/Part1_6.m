clear all
close all
clc;

load Data/PCAPCR.mat
f = FunctionCollection();
%% Part (a)
[Ux, Sx, Vx] = svd(X);
[Uxn, Sxn, Vxn] = svd(Xnoise);


% Plotting singular values of X and Xnoise
figure(1); 
subplot(1, 2, 1); 
hold on;
stem(1:10,diag(Sx),LineWidth=1.3); 
stem(1:10,diag(Sxn),LineWidth=1.0); 
grid on
set(gca,'fontsize', 14);
xlabel('Value Index'); 
ylabel('Singular Value magnitude');
title('Singular Values of X and Xnoise'); 
legend('X', 'Xnoise'); 
legend('$X$',' $X_{noise}$', 'Interpreter', 'latex');
hold off;

% Square error between SVs
subplot(1, 2, 2);
stem(1:size(Sx,2),(diag(Sx) - diag(Sxn)).^2,LineWidth=1.0,Color='g'); 
grid on
set(gca,'fontsize', 14);
xlabel('Value Index'); 
ylabel('Square Error');
title('Square Error between X and Xnoise');


%% Part (b)

[U, S, V] = svds(Xnoise,3);
Xdenoised = U*S*V';

% Square error between SVs
figure;
stem(1:10,f.MSE(X,Xnoise),LineWidth=1.5,Color='b'); 
hold on
stem(1:10,f.MSE(X,Xdenoised),LineWidth=1.0,Color='r'); 
grid on
set(gca,'fontsize', 14);
xlabel('Value Index'); 
ylabel('Mean Squared Error');
title('MSE between X and Xnoise (corrupted and denoised)');
legend('$X$, $X_{noise}$', '$X$, $\tilde{X}_{noise}$', 'Interpreter', 'latex');


%% Part (c)

B_ols = (Xnoise' * Xnoise) \ (X' * Y);
B_pcr = bPCR(Y, Xnoise, 3);

% Estimation and Test results
Y_ols = Xnoise * B_ols;
Y_pcr = Xnoise * B_pcr;
Y_test_ols = Xtest * B_ols;
Y_test_pcr = Xtest * B_pcr;

% Plotting squared errors
figure;
subplot(1, 2, 1);
hold on;
stem(1:5,f.MSE(Y,Y_ols),LineWidth=2);
stem(1:5,f.MSE(Y,Y_pcr),LineWidth=1.5);
xlabel('Value Index');
ylabel('Mean Squared Error');
title('MSE for Y');
legend('OLS', 'PCR');
grid minor;
hold off;

subplot(1, 2, 2);
hold on;
stem(1:5,f.MSE(Ytest,Y_test_ols),LineWidth=2);
stem(1:5,f.MSE(Ytest,Y_test_pcr),LineWidth=1.5);
xlabel('Value Index');
ylabel('Mean Squared Error');
title('MSE for Ytest ');
legend('OLS', 'PCR');
grid minor;
hold off;

%% Part (d)

ols_error = [];
pcr_error = [];

for i = 1:1000
    [Yhat_ols,Y_ols] = regval(B_ols);
    [Yhat_pcr,Y_pcr] = regval(B_pcr);
    ols_error = [ols_error; MSE(Y_ols,Yhat_ols)]; 
    pcr_error = [pcr_error; MSE(Y_pcr,Yhat_pcr)]; 
end


ols_error = mean(ols_error, 1);
pcr_error = mean(pcr_error, 1);

figure;
hold on;
stem(1:5,ols_error,LineWidth=2);
stem(1:5,pcr_error,LineWidth=1.5); 
xlabel('Value Index');
ylabel('MSE')
title('Ensemble Averaged Square Error Between Estimate & Original Data');
legend('OLS', 'PCR');
grid minor;
hold off;

%%
function B = bPCR(y,x,k)
    [U, S, V] = svds(x, k); % Retain only k largest components
    S_inv = diag(1 ./ diag(S)); % Compute pseudoinverse
    B = (V * S_inv * U') * y;
end
function error = MSE(x1,x2)
    %Mean Squared Error
    error = mean((x1 - x2).^2,1);
end
