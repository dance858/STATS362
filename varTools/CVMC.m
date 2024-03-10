function [losses, valueAtRiskEstimate, lb, ub, standardDeviation] = CVMC(m, nSamples, dt, Sigma_S, ...
    portfolioTheta, deltas, gammas,  lossFunc, threshold)
% m = number of market prices
% nSamples = number of replications to use for the estimation
% Sigma_S = covariance matrix of market moves
% portfolioTheta = portfolio theta
% deltas = portfolio deltas
% gammas = portfolio gammas, diagonal elements (off-diagonal elements are 
%          assumed to be 0).
% lossFunc = function that evaluates the portfolio loss given market moves dS.
% threshold = VAR threshold

CTilde = chol(Sigma_S, 'lower'); % CTilde*CTilde' = Sigma_S
assert(norm(CTilde*CTilde' - Sigma_S) < 1e-8);
[U, lambdas] = eig(-0.5*CTilde'*diag(gammas)*CTilde);
assert(norm(U*lambdas*U' + 0.5*CTilde'*diag(gammas)*CTilde) < 1e-8);
lambdas = diag(lambdas);
C = CTilde*U;
b = -C'*deltas;
a = -portfolioTheta*dt; 

% Generate Z
Z = randn(m, nSamples);
Q = a + b'*Z + lambdas'*(Z.^2);
dS = C*Z;
assert(isreal(Q));

losses = zeros(nSamples, 1);
for i = 1:nSamples
    losses(i, 1) = lossFunc(dS(:, i));
end

fX = (losses' > threshold);
hX = (Q > threshold);
assert(isreal(fX));
assert(isreal(hX))
assert(~isnan(sum(fX)))
assert(~isnan(sum(hX)))

% Prob(Q <= threshold)
myIntegrand = @(u) integrand(u, threshold, 1000, a, b, lambdas);
prob = (1/pi) * (integral(myIntegrand, 0, 40));
assert(isreal(prob));
assert(~isnan(prob))

beta = cov(fX, hX)/std(hX)^2;
beta = beta(1, 2);
valueAtRiskEstimate = mean(fX) - beta*(mean(hX) - (1 - prob));
standardDeviation = std(fX - beta*hX);
lb = valueAtRiskEstimate - 2.58*standardDeviation/sqrt(nSamples);
ub = valueAtRiskEstimate + 2.58*standardDeviation/sqrt(nSamples);

end


