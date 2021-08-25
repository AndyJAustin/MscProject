function p_Nmin1 = HermiteInterpolation(xk, fk, x, alpxk, alpx)
% Implements barycentric Hermite interpolation like in section 3.3 of Weideman and Reddy's
% MATLAB differentiation matrix suite. ACM TOMS, Vol. 26, pp. 465--519
% (2000). 
%f is of size n
%x is the point to be evaluated

%n = max(size(f));
%N = 2*n + 1;

%x_j = herroots(N);
%Scale the grid and the weights if the extremal values are >10^-6

%alpha = @(x) exp(-x.^2/2);
N = length(xk);
H_prime = @(x) 2*N*HermitePolynomial(N-1, x); %factor of 2 for physicists

%Setting up u_0 from f data (using spatial symmetry and knowing a0)
%u_0 = zeros(N,1);
%u_0(n+1) = a0;
%u_0(n+2:N) = f;
%u_0(1:n) = flipud(f);


phi = cell(N,1);
for j=1:N
    phi{j} = @(x) HermitePolynomial(N, x)./(H_prime(xk(j)).*(x - xk(j)));
end %phi{j}(x)=0 if x = x_j = 0. Do not let x coincide with the Hermite root at zero

p_Nmin1 = 0;

for j=1:N
    p_Nmin1 = p_Nmin1 + (alpx./alpxk(j)).*phi{j}(x).*fk(j);
end

end

