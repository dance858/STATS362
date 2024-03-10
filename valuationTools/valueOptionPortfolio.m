function [portfolioPrice, portfolioDeltas, portfolioGammas] = ...
    priceOptionPortfolio(St, dS, t, dt, T, r, KCall, KPut, volatility)
% St = column vector of m market prices at time t
% dt = risk-measurement horizon
% dS = column vector of change in S over interval dt
% T = expiration (we assume all stocks have the same expiration date).
% KCall = m x nCall, strike prices for call options
% KPut = m x nPut, strike prices for put options
% volatility = we assume all underlying assets have the same volatility

nAssets = length(St);
nCall = size(KCall, 2);
nPut = size(KPut, 2);
assert(nAssets == size(KCall, 1));
assert(nAssets == size(KPut, 1));

portfolioPrice = 0;
portfolioDeltas = zeros(nAssets, 1); % 
portfolioGammas = zeros(nAssets, 1); % the diagonal Gammas

for i = 1:nAssets
    % compute the price of all call options on asset i at time t + dt
    % given the price St at time t and the change dS over the time interval
    % dt.
    for j = 1:nCall
       [C, delta, gamma] = blackScholesCall(St(i) + dS(i), T, t + dt, r, KCall(i, j), volatility);
       portfolioPrice = portfolioPrice + C;
       portfolioDeltas(i) = portfolioDeltas(i) + delta;
       portfolioGammas(i) = portfolioGammas(i) + gamma;
    end

    % compute the price of all put options on asset i at time t + dt
    % given the price St at time t and the change dS over the time interval
    % dt.
    for j = 1:nPut
       [P, delta, gamma] = blackScholesPut(St(i) + dS(i), T, t + dt, r, KPut(i, j), volatility);
       portfolioPrice = portfolioPrice + P;
       portfolioDeltas(i) = portfolioDeltas(i) + delta;
       portfolioGammas(i) = portfolioGammas(i) + gamma;
    end
end
end