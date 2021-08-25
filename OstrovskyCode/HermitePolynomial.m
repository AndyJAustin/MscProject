function [ h ] = HermitePolynomial (n, x)
% Function that evaulates nth order Hermite polynomial H_n(x)
% Input
% -----
% * n: Order of polynomial
% * x: Input data (scalar or array) at which to evaluate polynomial.
%      - If no x input; function will return an array of the coefficients
%        of the polynomial.
%
% Output
% ------
% * h: Data of Hermite polynomial at x points, or the coefficients.

h = zeros(1,n+1);   %Set up arrays
n_fact = factorial(n);
m = 0:floor(n/2);

%Returning coeffficients.
%Phycisist's polynomials:
h(2*m+1) = n_fact .* (-1).^m ./ (factorial(m) .* factorial(n-2.*m)) .* 2.^(n-2.*m);
%Probabilist's polynomials:
%h(2*m+1) = 2^(-.5*n).* n_fact .* (-1).^m ./ (factorial(m) .* factorial(n-2.*m)) .* 2.^(.5.*(n-2.*m));

if exist('x', 'var') % If x is an input we evaluate H-polynomial at these points
    h = polyval(h, x);
end
end 