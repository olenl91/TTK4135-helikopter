function [c, ceq] = mycon_alt(x)
% system parameters
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;
N = 40;

lambda = x(1 : 6 : 6*N);
e = x(5 : 6 : 6*N);

c = zeros(N, 1);

for n = 1:N
    c(n) = alpha * exp(-beta .* (lambda(n) - lambda_t).^2) - e(n);
end

ceq = [];
%mycon = @(x) [c(x), []];

end

