function [u_lj, u_lj_t, G, dGdT, dGdc, momentum, tdata] = FIOstrovskyUneven(u0, T, c, M, J, x)
% Function that integrates he Ostrovsky equation forward in time in order
% to compute the objective function and its derivatives with respect to c
% and T. Returns a set of u0 values that are needed for the adjoint
% integration, among other useful data for analysis and visualisation.
% Includes a dampening region to prevent dispersion reentering the periodic
% integration domain. 
%
% Uses alternative formulation of objective that relies on integration over
% an entire period.
%
% Uses a Fourier pseudo-spectral integrating factor method similar to that
% presented in L. Trefethen, "Spectral Methods in MATLAB", SIAM,
% Philadelphia, 2000. 
%
% Input
% -----
% * u0: initial wavefunction on L Fourier nodes
% * T: Half-period of wave.
% * c: Groupspeed of wave.
% * M: Number of time steps, must be large enough for accuracy and
%      stability. Cost is linear in M.
% * J: Number of slices that the time steps are split into to return a
%      reduced matrix of u0 "milepost" values for adjoint integration, where the
%      values of u0(x,t) in between milepost times are interpolated when
%      needed in adjoint RK4.
% * x: Evenly spaced Fourier nodes associated with u0 data.
%
% Output
% ------
% * u_lj: Matrix containing Milepost values of u(x,t), named such as there
%         are J+1 mileposts denoted by indices j. Each column is a
%         wavefunction u0(x,t_j) at a milepost time, and the rows denote
%         the Fourier nodes with index l.
% * u_lj_t: Matrix containing milepost values of the temporal derivative
%         of the solution u.
% * G: Computed value of the objective function.
% * dGdT: Derivative of the objective function wrt half-period.
% * dGdc: Derivative of the objective function wrt groupspeed.
% * momentum: Computed value of momentum flux p = 1/2 integral[u^2]dx
% * tdata: Array of values of time at each time step.
%
% The degree of rotation may be controlled with the variable epsilon. Kept
% inside the function not to clutter the function calls.

dt = 2*T/M;
Q = M/J; % number of nodes per slice
L = length(u0);

x_min = x(1);
x_max = - x_min;
xfac = (x_max - x_min)/(2*pi);

u = u0;
p_integrand = [u; u(1)].^2;
momentum = .5*trapz([x;x_max], p_integrand);

% Wavenumbers
kap = [0:L/2-1 0 -L/2+1:-1]'; %zero in middle to get rid of Nyquist frequency
k = kap/xfac;
epsilon = 1;
iom = 1i*k*c + 1i*k.^3 - epsilon^2*1i./(k+eps);

%Damping
kapwid = 1;
kap0 = x_min + 3*kapwid;
kapamp = 0.05;
kappa = kapamp*exp(-((x - kap0)/kapwid).^2);

%Initializing variables used in calculating derivatives of G_tot
j = 1;
u_lj = zeros(L, J + 1);
u_lj_t = zeros(L,J+1);
u_lj(:,j) = u0;

% Solve PDE over T
v = fft(u);
E = exp(dt*iom/2);
E2 = E.^2;
g = -0.5i*dt*k;
tdata = zeros(J+1,1);
tdata(1) = 0;
for m = 1:M
   t = m*dt;
   u = real( ifft(     v    ) );
   a = g.*fft(u.^2)- fft(kappa .* u);
   u = real( ifft(E.*(v + a/2)) );
   b = g.*fft(u.^2)- fft(kappa .* u);
   u = real( ifft(E.*v + b/2));
   c = g.*fft(u.^2)- fft(kappa .* u);
   u = real( ifft(E2.*v + E.*c) );
   d = g.*fft(u.^2)- fft(kappa .* u);
   v = E2.*v + (E2.*a + 2*E.*(b+c) + d)/6;
   v(L/2+1) = 0;
   v(L/2+2) = 0;
   %save components for du/dt at nodes for use in adjoint

   if mod(m, round(Q)) == 0 && m ~= 1
        % add values to u_lj
        j = j + 1;
        u_lj(:, j) = u;

        u_lj_t(:,j) = real(ifft(E2.*a/dt));
        tdata(j) = t; %for visualisation
        if m == M
            u_l2T_x = real(ifft(1i*k.*v));
        end
   
   elseif m == 1
       u_lj_t(:,j) = real(ifft(E2.*a/dt));
       
   end
    
end

% calculate G and its derivatives wrt T and c
u_l2T = u_lj(:,J+1);
u_l2T_t = u_lj_t(:,J+1);
% -- G --
G_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).^2;
G = .5*trapz([x;x_max], G_integrand);

% -- dG/dT --
dGdT_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).*[u_l2T_t; u_l2T_t(1)];
dGdT = 2*trapz([x;x_max], dGdT_integrand);

% -- dG/dc --
dGdc_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).*[u_l2T_x; u_l2T_x(1)];
dGdc =  -T*2* trapz([x;x_max], dGdc_integrand);
end
