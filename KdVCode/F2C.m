function [dGduC] = F2C(dGduF, Symm)
% Function that converts the derivative of G with respect to u0 on evenly
% spaced Fourier nodes to the derivative of G with respec to Fourier mode
% coefficients by means of Fourier transformation of the individual
% componenents.
%
% Input
% -----
% * dGduF: Derivative of G with respect to u0 on L evenly spaced Fourier
%          nodes 
% * Symm: Logical value for symmetry. 
%
% Output
% ------
% * dGduC: Derivative of G with respect to u0 on L/2-1 or L/4 Fouier
%          coefficinets c0, a1, a2,..., {b1, b2,...}, where c_k = a_k + i
%          b_k.

L = length(dGduF);
dGdv = fft(dGduF);
if isequal(Symm, true)
    dGduC = zeros(L/4,1); 
    dGduC(1) = 2*pi*dGdv(1);
    dGduC(2:L/4) = 4*pi*real(dGdv(2:L/4));
else
    dGduC = zeros(L/2-1,1); 
    dGduC(1) = 2*pi*dGdv(1);
    dGduC(2:L/4) = 4*pi*real(dGdv(2:L/4));
    dGduC(L/4+1:L/2-1) = 4*pi*imag(dGdv(2:L/4));

end

end