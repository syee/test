clear;
clc;
rng(1111191);

%*************************************************************************%
[beta_hat,sigma_beta_hat,sigma_beta_tilde]=simulateFEest(5,1);

%1a
figure(1)
hist(beta_hat); 
xlabel('Beta estimates')
ylabel('Count')
title('Beta estimate values (T=5)')
saveas(figure,'figure1.jpeg')

figure(2)
hist(sigma_beta_hat);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, robust to autocorrelation (T=5)')

figure(3)
hist(sigma_beta_tilde);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, assume no autocorrelation (T=5)')

%1b
true_sd = std(beta_hat);

%1c
bias_sd_hat = true_sd - mean(sigma_beta_hat);
bias_sd_tilde = true_sd - mean(sigma_beta_tilde);

std_sd_hat = std(sigma_beta_hat);
std_sd_tilde = std(sigma_beta_tilde);

rmse_sd_hat_5 = sqrt(bias_sd_hat^2 + std_sd_hat^2);
rmse_sd_tilde_5 = sqrt(bias_sd_tilde^2 + std_sd_tilde^2);


%*************************************************************************%
%1d
%Repeat the simulation with T = 10.
[beta_hat_d,sigma_beta_hat_d,sigma_beta_tilde_d]=simulateFEest(10,1);

figure(4)
hist(beta_hat_d); 
xlabel('Beta estimates')
ylabel('Count')
title('Beta estimate values (T=10)')

figure(5)
hist(sigma_beta_hat_d);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, robust to autocorrelation (T=10)')

figure(6)
hist(sigma_beta_tilde_d);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, assume no autocorrelation (T=10)')

true_sd_d = std(beta_hat_d);

bias_sd_hat_d = true_sd_d - mean(sigma_beta_hat_d);
bias_sd_tilde_d = true_sd_d - mean(sigma_beta_tilde_d);

std_sd_hat_d = std(sigma_beta_hat_d);
std_sd_tilde_d = std(sigma_beta_tilde_d);

rmse_sd_hat_10 = sqrt(bias_sd_hat_d^2 + std_sd_hat_d^2);
rmse_sd_tilde_10 = sqrt(bias_sd_tilde_d^2 + std_sd_tilde_d^2);


%*************************************************************************%
%Repeat the simulation with T = 20.
[beta_hat_d2,sigma_beta_hat_d2,sigma_beta_tilde_d2]=simulateFEest(20,1);

figure(7)
hist(beta_hat_d2); 
xlabel('Beta estimates')
ylabel('Count')
title('Beta estimate values (T=20)')

figure(8)
hist(sigma_beta_hat_d2);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, robust to autocorrelation (T=20)')

figure(9)
hist(sigma_beta_tilde_d2);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, assume no autocorrelation (T=20)')

true_sd_d2 = std(beta_hat_d2);

bias_sd_hat_d2 = true_sd_d2 - mean(sigma_beta_hat_d2);
bias_sd_tilde_d2 = true_sd_d2 - mean(sigma_beta_tilde_d2);

std_sd_hat_d2 = std(sigma_beta_hat_d2);
std_sd_tilde_d2 = std(sigma_beta_tilde_d2);

rmse_sd_hat_20 = sqrt(bias_sd_hat_d2^2 + std_sd_hat_d2^2);
rmse_sd_tilde_20 = sqrt(bias_sd_tilde_d2^2 + std_sd_tilde_d2^2);


%*************************************************************************%
%1f
% Change the value of beta to 314
[beta_hat_f,sigma_beta_hat_f,sigma_beta_tilde_f]=simulateFEest(5,314.15926);

true_sd_f = std(beta_hat_f);

figure(10)
hist(beta_hat_f); 
xlabel('Beta estimates')
ylabel('Count')
title('Beta estimate values (T=5), higher true beta value')

figure(11)
hist(sigma_beta_hat_f);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, robust to autocorrelation (T=5)')

figure(12)
hist(sigma_beta_tilde_f);
xlabel('Standard error of beta values')
ylabel('Count')
title('Standard error of beta, assume no autocorrelation (T=5)')

bias_sd_hat_b = true_sd_f - mean(sigma_beta_hat_f);
bias_sd_tilde_f = true_sd_f - mean(sigma_beta_tilde_f);

std_sd_hat_b = std(sigma_beta_hat_f);
std_sd_tilde_f = std(sigma_beta_tilde_f);

rmse_sd_hat_b = sqrt(bias_sd_hat_f^2 + std_sd_hat_f^2);
rmse_sd_tilde_b = sqrt(bias_sd_tilde_f^2 + std_sd_tilde_f^2);

function [beta_hat,sigma_beta_hat,sigma_beta_tilde] = simulateFEest(t,beta_val)
    N = 500;
    T = t;
    nreplic = 1000;
    beta = beta_val;

    beta_hat = zeros(nreplic,1);

    %Allows for possible serial correlation
    sigma_beta_hat =  zeros(nreplic,1);

    %Assuming that serial correlation is 0
    sigma_beta_tilde =  zeros(nreplic,1);

    for replic = 1:1:nreplic

        %Create a micro panel from the standard normal.
        X = randn(T,N);
        %Sample covariance 
        u = abs(X).*randn(T,N);
        y = X + u;

        %Demeaned variables
        ydm = y - repmat(mean(y),T,1);
        udm = u - repmat(mean(u),T,1);
        Xdm = X - repmat(mean(X),T,1);

        beta_hat(replic) = beta + (sum(sum(Xdm.^2,2),1)).^(-1) * sum(sum(Xdm.*udm,2),1);
 
        % residual from fixed effect regression
        e_hat = ydm - Xdm.*beta_hat(replic); 

        S_xx = sum(Xdm.^2,'all');
        
        %Use formula 3 to calculate the standard error of the estimator.
        sigma_beta_hat(replic) = sqrt(S_xx^(-2) * sum(sum(Xdm.*e_hat,1).^2,2));

        %Use formula 4 to calculate the standard error of the estimator.
        sigma_beta_tilde(replic) = sqrt(S_xx^(-2) * sum(sum(Xdm.^2.*e_hat.^2,1),2));
    end

end
