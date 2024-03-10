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
deltas = -deltas;
gammas = -gammas;
theta = -theta;

nSamples = 50000;
all_thresholds = [80, 100, 120, 140, 160, 180, 200];
lbs_plainMC = zeros(length(all_thresholds), 1);
ubs_plainMC = zeros(length(all_thresholds), 1);
VaR_plainMC = zeros(length(all_thresholds), 1);
stds_plainMC = zeros(length(all_thresholds), 1);
lbs_IS = zeros(length(all_thresholds), 1);
ubs_IS = zeros(length(all_thresholds), 1);
VaR_IS = zeros(length(all_thresholds), 1);
stds_IS = zeros(length(all_thresholds), 1);
lbs_CV = zeros(length(all_thresholds), 1);
ubs_CV = zeros(length(all_thresholds), 1);
VaR_CV = zeros(length(all_thresholds), 1);
stds_CV = zeros(length(all_thresholds), 1);

for i = 1:length(all_thresholds)
    threshold = all_thresholds(i);
    fprintf("Simulating threshold = %i. \n", threshold);

    [lossesPlain, VARPlainEst, lbPlain, ubPlain, stdPlain] = ...
        plainMC(m, nSamples, Sigma_S, lossFunc, threshold);

    lbs_plainMC(i) = lbPlain;
    ubs_plainMC(i) = ubPlain;
    VaR_plainMC(i) = VARPlainEst;
    stds_plainMC(i) = stdPlain;
    
    [lossesIS, VaRISEst, lbIS, ubIS, stdIS] = ...
        ISMC(m, nSamples, dt, Sigma_S, theta, deltas, gammas,  lossFunc, threshold);
    
    lbs_IS(i) = lbIS;
    ubs_IS(i) = ubIS;
    VaR_IS(i) = VaRISEst;
    stds_IS(i) = stdIS;

    [lossesCV, VaRCVEst, lbCV, ubCV, stdCV] = ...
        CVMC(m, nSamples, dt, Sigma_S, theta, deltas, gammas,  lossFunc, threshold);
    
    lbs_CV(i) = lbCV;
    ubs_CV(i) = ubCV;
    VaR_CV(i) = VaRCVEst;
    stds_CV(i) = stdCV;
end 

if saveData
   save('VaRVsThresholdData.mat')
end

%%
%all_thresholds = all_thresholds(2:end);
%VaR_plainMC = VaR_plainMC(2:end);
%lbs_plainMC = lbs_plainMC(2:end);
%ubs_plainMC = ubs_plainMC(2:end);
%VaR_IS = VaR_IS(2:end);
%lbs_IS = lbs_IS(2:end);
%ubs_IS = ubs_IS(2:end);
%VaR_CV = VaR_CV(2:end);
%lbs_CV = lbs_CV(2:end);
%ubs_CV = ubs_CV(2:end);



hold on; 
line1 = plot(all_thresholds, VaR_plainMC, 'k-o', 'LineWidth', 2);
line2 = plot(all_thresholds, lbs_plainMC, 'k--', all_thresholds, ubs_plainMC, 'k--', 'LineWidth', 1.5); 
line3 = plot(all_thresholds, VaR_IS, 'b-x', 'LineWidth', 2);
line4 = plot(all_thresholds, lbs_IS, 'b--', all_thresholds, ubs_IS, 'b--', 'LineWidth', 1.5); 
%line5 = plot(all_thresholds, VaR_CV, 'r-', 'LineWidth', 2);
%line6 = plot(all_thresholds, lbs_CV, 'r:', all_thresholds, ubs_CV, 'r:', 'LineWidth', 1.5); 

% Adding labels and legend
xlabel('Thresholds');
ylabel('VaR')
legend('Plain estimate', 'Plain conf. interval', 'IS estimate', 'IS conf. interval')%, ...
    %'CV estimate')%, 'CV LB', 'CV UB', 'Location', 'best');
%legend('Plain estimate', 'Plain LB', 'Plain UB', 'IS estimate', 'IS LB', ...
%        'IS UB', 'CV estimate', 'CV LB', 'CV UB', 'Location', 'best');

% Display grid
grid on;

% Release the hold
hold off;


