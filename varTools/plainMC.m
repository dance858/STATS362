function [losses, valuteAtRiskEstimate, lb, ub, stdPlain] = plainMC(m, nSamples, Sigma_S, lossFunc, threshold)
% m = number of market prices
% nSamples = number of replications to use for the estimation
% Sigma_S = covariance matrix of market moves
% lossFunc = function that evaluates the portfolio loss given market moves dS.
% threshold = VAR threshold

% generate samples 
L = chol(Sigma_S, 'lower'); % Sigma_S = L*L'
assert(norm(Sigma_S - L*L') < 1e-8);
dS = L*randn(m, nSamples);

losses = zeros(nSamples, 1);
for i = 1:nSamples
    losses(i, 1) = lossFunc(dS(:, i));
end

fX = (losses > threshold);
valuteAtRiskEstimate = mean(fX);
stdPlain = std(fX);

lb = valuteAtRiskEstimate - 2.58*stdPlain/sqrt(nSamples);
ub = valuteAtRiskEstimate + 2.58*stdPlain/sqrt(nSamples);
end