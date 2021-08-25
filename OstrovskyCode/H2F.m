function [u0, x, c, T] = H2F(q, b, L, Symm)
% Function that converts the q into evenly spaced Fourier nodes, and the
% individual values of c and T using the information about the
% representation of q and the number of Fourier points to be interpolated
% onto. Padds each side of the Hermite points to zero so is only suitable
% with wavepackets and an appropriate scaing factor b, such that the values
% at the extremal nodes are negligable. Makes use of piecewise cubic
% Hermite interpolation.
%
% Input
% -----
% * q: The desccent coordinates vector (u0,c,t)
% * b: Scaling parameter x = b x_tilde.
% * L: Number of Fourier nodes.
% * Symm: Logical value for symmetry. 
%
% Output
% ------
% * u0: Nodal values of u0 at L Fourier points.
% * x: x-values corresponding to u0 data.
% * c: Groupspeed of wavepacket.
% * T: Half-period of wavepacket.

d = length(q);
c = q(d-1);
T = q(d);

% ----- Set Up New u0 Array -----

if isequal(Symm, true)
    n = d - 3;
    N = 2*n + 1;
    u0_tot = zeros(N, 1);
    u0_tot(1:n+1) = q(1:n+1); u0_tot(n+2:N) = flipud(q(1:n));
else
    N = d - 2;
    u0_tot = q(1:N);
end

% ----- Compute u0 and x -----

x_j = herroots(N);
x_scaled = x_j./b;
x_min = 1.05*x_scaled(1);
x_max = -x_min;
xfac = (x_max - x_min)/(2*pi);
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = xfac*xi;
u0 = interp1([x_min;x_scaled;x_max],[0;u0_tot;0],x,'pchip');

end