function [dGduH] = F2H(dGduF, b, N, Symm)
% Function that converts the derivative of G with respect to u0 on evenly
% spaced Fourier nodes to the derivative of G with respect to Hermite nodes
% by means of Fourier interpolant. Takes care of symmetry or lack thereof.
%
% The Fourier interpolant (fourint.m) used is taken from Weideman and
% Reddy's MATLAB differentiation matrix suite. ACM TOMS, Vol. 26, pp.
% 465--519 (2000). https://appliedmaths.sun.ac.za/~weideman/research/differ.html
%
% Input
% -----
% * dGduF: Derivative of G with respect to u0 on L evenly spaced Fourier
%          nodes 
% * b: Scaling factor x = b x_tilde.
% * N: Number of Hermite nodes to be interpolated back to.
% * Symm: Logical value for symmetry. 
%
% Output
% ------
% * dGduH: Derivative of G with respect to u0 on N or (N-1)/2 Hermite
%          nodes. 

x_j = herroots(N);
x_scaled = x_j./b;
x_min = 1.05*x_scaled(1); % 2*x_scaled(1);
x_max = -x_min;      xfac = (x_max - x_min)/(2*pi);
%L = length(dGduF);
%dGduF(1:L/4-1) = 0; %Only do this when we use a much larger x-domain for
%dGduF(L/4+1:L) = 0; %integration and can ignore effeccts outside the
                     %envelope region
w = fourint(dGduF, x_j/(b*xfac));
if isequal(Symm, true)
    n = (N-1)/2;
    dGduH = w(1:n+1);
else
    dGduH = w;
end
end