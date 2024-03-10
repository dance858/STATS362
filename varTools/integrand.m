function result = integrand(u, x, y, a, b, lambdas)

theta = 1i*u;
psi = a.*theta + 0.5*sum((theta.^2 .* b.^2)./(1-2.*theta.*lambdas) - log(1 - 2.*theta.*lambdas), 1);
%assert(isreal(psi));
phiHat = exp(psi);
result = real(phiHat .* ((exp(1i.*u*y) - 1)./(1i.*u)) .* exp(-1i.*u*x));
end