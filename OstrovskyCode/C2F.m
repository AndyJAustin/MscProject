function [u0, x] = C2F(v0, xfac, Symm)
% Function that converts the Fourier mode coefficients into evenly spaced
% Fourier nodes by means of an inverse Fourier transform.
%
% Input
% -----
% * v0: The Fourier coefficients arranged as c0, a1, a2,..., {b1, b2,...}
%       of length L/2-1 or L/4 depending on whether imaginary coefficints
%       are excluded.
% * xfac: Scaling parameter xfac=(xmax-xmin)/2pi.
% * Symm: Logical value for symmetry. 
%
% Output
% ------
% * u0: Nodal values of u0 at L Fourier points.
% * x: x-values corresponding to u0 data.

f = length(v0);
if isequal(Symm, true)
    L = f*4;  %f = L/4
    v0 = zeros(L, 1);
    v0(1:L/4) = q(1:L/4);
    v0(3*L/4+2:L) = flipud(v0(2:L/4));
    u0 = real(ifft(v0));
else
    L = (f+1)*2; %f = L/2-1
    v0 = zeros(L,1);
    v0(1:L/4) = q(1:L/4);
    v0(2:L/4) = v0(2:L/4) + 1i*q(L/4+1:f); 
    v0(3*L/4+2:L) = flipud(conj(v0(2:L/4)));
    u0 = real(ifft(v0));
    
end
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = xfac*xi;

end