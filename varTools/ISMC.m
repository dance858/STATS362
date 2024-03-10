function [losses, valueAtRiskEstimate, lb, ub, standardDeviation] = ISMC(m, nSamples, dt, Sigma_S, ...
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

% Pick value of theta
gFunction = @(theta) a + sum(theta*(b.^2).*(1-theta*lambdas)./((1-2*theta*lambdas).^2) + lambdas./(1-2*theta*lambdas)) - threshold; 
tic;
theta = fzero(gFunction, 0);
toc;
phiTheta = a*theta + 0.5*sum((theta^2*b.^2)./(1-2*theta*lambdas) - log(1 - 2*theta*lambdas));
assert(isreal(phiTheta));
assert(isreal(theta));
assert(~isnan(theta))
assert(max(theta*lambdas) < 0.5)

% Generate Z
muTheta = theta*b./(1-2*lambdas*theta);
sigmasSquaredTheta = 1./(1-2*theta*lambdas); % Should never be negative...
Z = muTheta + diag(sqrt(abs(sigmasSquaredTheta)))*randn(m, nSamples);
assert(isreal(Z));
% Evaluate Q
Q = a + b'*Z + lambdas'*(Z.^2);
dS = C*Z;
assert(isreal(Q));

losses = zeros(nSamples, 1);
for i = 1:nSamples
    losses(i, 1) = lossFunc(dS(:, i));
end

fXIs = exp(-theta*Q + phiTheta).*(losses' > threshold);

standardDeviation = std(fXIs);
valueAtRiskEstimate = mean(fXIs);
lb = valueAtRiskEstimate - 2.58*standardDeviation/sqrt(nSamples);
ub = valueAtRiskEstimate + 2.58*standardDeviation/sqrt(nSamples);
end