function [C, delta, gamma, theta] = blackScholesCall(St, T, t, r, K, volatility)
% Black-Scholes for call option
% St = price of underlying asset at time t
% T = expiration date
% t = time 
% r = annualized risk-free interest rate, continuously compunded
% K = strike price
% volatility = standard deviation of stock returns

dPlus = 1/(volatility*sqrt(T-t))*(log(St/K) + (r + volatility^2/2)*(T-t));
dMinus = dPlus - volatility*sqrt(T-t);

temp = normcdf(dMinus);
delta = normcdf(dPlus);
gamma = 1/sqrt(2*pi)*exp(-dPlus^2/2)/(St*volatility*sqrt(T-t));
theta = -St*normpdf(dPlus)*volatility/(2*sqrt(T-t)) - r*K*exp(-r*(T-t))*temp;
C = delta*St - temp*K*exp(-r*(T-t));
end