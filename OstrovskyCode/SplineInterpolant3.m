function [u0, x] = SplineInterpolant3(u_0, b, L)
%Used in FIOstrovskyDiffMat to calculate momentum, objective and its
%derivatives wrt c and T. 
%
% We do not proceed with differentiation matrix approach, and instead opt
% for an integrating factor method, as the computational cost outweigh the
% marginal benefits offered. 

N = length(u_0);
x_j = herroots(N);
x_scaled = x_j./b;
x_min = 1.05*x_scaled(1); %Extra to make sure the x values extend further than extremal H nodes
x_max = -x_min;
xfac = (x_max - x_min)/(2*pi);
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = 0.5*(x_max + x_min) + xfac*xi;
u0 = interp1([x_min;x_scaled;x_max],[0;u_0;0],x,'pchip');
%Padded with zeros at end to stop edge effects
end