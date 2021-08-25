function [dGdu, sdata, wdata] = AIOstrovskyEvenxperiodic(u_lj, u_lj_t, T, M, J, x)
% Function that performs an adjoint integration to compute the varational
% derivative of the objective function with respect to the initial
% condition u0. Uses information from the forward integration during the
% to comopute values of u(x,t) inbetween time steps through cubic Hermite
% Interpolation and implements a Fourier pseudo-spectral integrating
% factor method similar to that presented in L. Trefethen, "Spectral
% Methods in MATLAB", SIAM, Philadelphia, 2000. 
%
% Used for a spatially periodic Ostrovsky solution approximation
%
% Input
% -----
% * u_lj: Matrix containing Milepost values of u(x,t), named such as there
%         are J+1 mileposts denoted by indices j. Each column is a
%         wavefunction u0(x,t_j) at a milepost time, and the rows denote
%         the Fourier nodes with index l.
% * u_lj_t: Matrix containing milepost values of the temporal derivative
%         of the solution u for interpolation.
% * T: Half-period of wave.
% * M: Number of time steps, must be large enough for accuracy and
%      stability. Cost is linear in M.
% * J: Number of mileposts minus 1.
% * x: Evenly spaced Fourier nodes associated with u0 data.
%
% Output
% ------
% * dGdu: Variational derivative at Fourier nodes.
% * sdata: "Reverse time" step data; s = T - t
% * wdata: Matrix of adjoint function w(x,s) at each milepost for
%          visualisation of integration.
%
% The degree of rotation may be controlled with the variable epsilon. Kept
% inside the function not to clutter the optimisation function calls.


%Set up initial parameters
L = size(u_lj, 1);
x_min = x(1);
x_max = -x_min;
xfac = (x_max - x_min)/(2*pi);
ds = T/M;
Q = M/J; % nodes per slice - must be integer

% Wavenumbers
kap = [0:L/2-1 0 -L/2+1:-1]';
k = kap/xfac;
epsilon = 1e-8;
iom = 1i*k.^3 - epsilon^2*1i./(k+eps);

%Set up initial values
u_high = u_lj(:, J + 1); % u(x,T)
u_low = u_lj(:, J);
u_t_high = u_lj_t(:, J + 1);
u_t_low = u_lj_t(:, J);
u_negx = zeros(L, 1);
u_negx(2:L) = flipud(u_high(2:L)); % u(-x,T)
u_negx(1) = u_high(1);
w_hat = fft(u_high - u_negx); %initial condition

% Cubic Hermite Interpolation:
% u(.,t_n + ph) = (1-p)u_n + pu_nplus1 - p(1-p){(1-2p)(u_nplus1-u_n) -
%(1-p)h du_n/dt + ph du_nplus1/dt}

h = T/J; % h is the timestep between mileposts
theta = 1;

%Solve adjoint equation over all discretised values of s
E = exp(-ds*iom/2);
E2 = E.^2;
g = 1i*ds*k;
sdata = zeros(1, J+1);
sdata(1) = 0;
wdata = 2*real(ifft(w_hat));
for m = 1:M
   s = m*ds;
   %Cubic Hermite interpolation
   u = (1-theta)*u_low+theta*u_high-theta*(1-theta)*((1-2*theta)*(u_high-u_low) ...
   - (1-theta)*h*u_t_low + theta*h*u_t_high);

    theta2 = theta - 1/(2*Q);
    
    u_half = (1-theta2)*u_low+theta2*u_high-theta2*(1-theta2)*((1-2*theta2)*(u_high-u_low) ...
   - (1-theta2)*h*u_t_low + theta2*h*u_t_high);

    theta = theta - 1/Q;
    
    u_one = (1-theta)*u_low+theta*u_high-theta*(1-theta)*((1-2*theta)*(u_high-u_low) ...
   - (1-theta)*h*u_t_low + theta*h*u_t_high);
    
    
   %s --> T-t
   w = real(ifft( w_hat ));
   a = g.*fft(u.*w);
     %s+ds/2 --> T-t-dt/2
   w = real( ifft(E.*(w_hat+a/2)) );
   b = g.*fft(u_half.*w);                        % 4th-order
     %s+ds/2 --> T-t-dt/2
   w = real( ifft(E.*w_hat + b/2) );
   c = g.*fft(u_half.*w);
     %s+ds --> T -t-dt
   w = real( ifft(E2.*w_hat+E.*c) );
   d = g.*fft(u_one.*w);
   w_hat = E2.*w_hat + (E2.*a + 2*E.*(b+c) + d)/6;
   %w_hat(1) = 0;
   w_hat(L/2 + 1) = 0;
   
   %Recalculate u_lm each time the adjoint integration passes a milepost
   if m < M && mod(m, round(Q)) == 0 
       theta = 1;
       j = J + 1 - round(m/Q); %index must be integer
       u_high = u_lj(:,j);
       u_low = u_lj(:,j-1);
       u_t_high = u_lj_t(:,j);
       u_t_low = u_lj(:,j-1);
       sdata(round(m/Q+1)) = s;
       wdata = [wdata 2*real(ifft(w_hat))];
   end
   if m == M
       wdata = [wdata 2*real(ifft(w_hat))];
       sdata(m/Q+1) = s;
   end
   
end

w = 2*real(ifft(w_hat));
%w((L/2+2:L)) = flipud(w(2:L/2)); 
w(2:L) = .5*(w(2:L) + flipud(w(2:L)));
dGdu = w;

end