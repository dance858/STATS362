% This script plots the VaR estimate and the confidence interval vs n.
clear; clc; addpath('valuationTools/', 'varTools/');

%% Problem data
rng(0)
m = 10; nCall = 10; nPuts = 6;
St = 100 + 50*(rand(m, 1) - 1); 
t = 0; T = 0.1; r = 0.05; 
KCall = St - 2*[-nCall/2:nCall/2];
KPut = St - 3*[-nPuts/2:nPuts/2];
dt = 0.04; volatility = 0.3;
Sigma_S = diag(volatility*sqrt(dt).*St).^2;
saveData = false;

%% Define loss function and threshold for VAR
[startPricePortfolio, deltas, gammas, theta] = ...
    priceOptionPortfolio(St, zeros(m, 1), t, 0, T, r, KCall, KPut, volatility);
lossFunc = @(dS) (priceOptionPortfolio(St, dS, t, dt, T, r, KCall, KPut, volatility) ...
                  - startPricePortfolio);
threshold = 120;
deltas = -deltas;
gammas = -gammas;
theta = -theta;


%% 
all_n = [10000:10000:20000];
lbs_plainMC = zeros(length(all_n), 1);
ubs_plainMC = zeros(length(all_n), 1);
VaR_plainMC = zeros(length(all_n), 1);
stds_plainMC = zeros(length(all_n), 1);
lbs_IS = zeros(length(all_n), 1);
ubs_IS = zeros(length(all_n), 1);
VaR_IS = zeros(length(all_n), 1);
stds_IS = zeros(length(all_n), 1);
lbs_CV = zeros(length(all_n), 1);
ubs_CV = zeros(length(all_n), 1);
VaR_CV = zeros(length(all_n), 1);
stds_CV = zeros(length(all_n), 1);


for i = 1:length(all_n)
    nSamples = all_n(i);
    fprintf("Simulating n = %i. \n", nSamples);

    [lossesPlain, VARPlainEst, lbPlain, ubPlain, stdPlain] = ...
        plainMC(m, nSamples, Sigma_S, lossFunc, threshold);

    fprintf("Plain MC done \n");
    lbs_plainMC(i) = lbPlain;
    ubs_plainMC(i) = ubPlain;
    VaR_plainMC(i) = VARPlainEst;
    stds_plainMC(i) = stdPlain;
    
    [lossesIS, VaRISEst, lbIS, ubIS, stdIS] = ...
        ISMC(m, nSamples, dt, Sigma_S, theta, deltas, gammas,  lossFunc, threshold);
    fprintf("Importance sampling done \n");
    
    lbs_IS(i) = lbIS;
    ubs_IS(i) = ubIS;
    VaR_IS(i) = VaRISEst;
    stds_IS(i) = stdIS;

    %% Control variate
    [lossesCV, VaRCVEst, lbCV, ubCV, stdCV] = ...
        CVMC(m, nSamples, dt, Sigma_S, theta, deltas, gammas,  lossFunc, threshold);
    
    fprintf("Control variate done \n");
    lbs_CV(i) = lbCV;
    ubs_CV(i) = ubCV;
    VaR_CV(i) = VaRCVEst;
    stds_CV(i) = stdCV;


end
if saveData
    save('estimateVsSamplesDataRealRun.mat')
end
%%
% Visualize confidence intervals
figure; hold on;

line1 = plot(all_n, VaR_plainMC, 'k-x', 'LineWidth', 2, 'MarkerSize', 8);
line2 = plot(all_n, lbs_plainMC, 'k--');
plot(all_n, ubs_plainMC, 'k--')
line3 = plot(all_n, VaR_IS, 'b-^', 'LineWidth', 2, 'MarkerSize', 8);
line4 = plot(all_n, lbs_IS, 'b--');
plot(all_n, ubs_IS, 'b--');
%line5 = plot(all_n, VaR_CV, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
%#line6 = plot(all_n, lbs_CV, 'r--');
%plot(all_n, ubs_CV, 'r--');




%fill([all_n, fliplr(all_n)], [lbs_plainMC', fliplr(ubs_plainMC')], ...
%    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%plot(all_n, VaR_IS, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%fill([all_n, fliplr(all_n)], [lbs_IS', fliplr(ubs_IS')], ...
%    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%plot(all_n, VaR_CV, '-o', 'LineWidth', 2, 'MarkerSize', 8);
%fill([all_n, fliplr(all_n)], [lbs_CV', fliplr(ubs_CV')], ...
%    'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

xlabel('Number of Samples');
ylabel('Estimate');
%legend([line1, line2, line3, line4, line5, line6], 'Plain estimate',  ...
%       'Plain conf. interval', 'IS estimate', 'IS conf. interval',  ...
%       'CV estimate', 'CV conf. interval', 'Location', 'Best');
legend([line1, line2, line3, line4], 'Plain estimate',  ...
       'Plain conf. interval', 'IS estimate', 'IS conf. interval', 'Location', 'Best');

grid on; 
hold off; 

%% Plot variance reduction vs samples

figure; hold on;
plot(all_n, (stds_plainMC./stds_IS).^2, '-o')
plot(all_n, (stds_plainMC./stds_CV).^2, '-x')
xlabel('Number of Samples');
ylabel('Improvement');
legend('IS', 'CV', 'Location', 'Best');
ylim([1, 12])
grid on


