%% Periodic Solutions
%  ------------------
%Set up periodic KdV IVs
L = 256;
r1 = 1.6;%1.1;
r2 = 2;% 1.5;       For T' to be small enough for stability, r2 and r3 must
r3 = 2.1; %4        be far enough apart. Check stability by integrating over long enough time
cg = (r1 + r2 + r3)/6; 
a = r2 - r1;
m = (r2 - r1)/(r3 - r1);
umin = .5*(r1 + r3 - r2);
K = ellipke(m);
k = ((r3 - r1)^(.5))/(4*3^.5*K);
chi = 1/(2*k); %Half of a period in x

x_max = chi;
x_min = -x_max;
xfac = (x_max - x_min)/(2*pi);
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = 0.5*(x_max + x_min) + xfac*xi;
%x = x + (pi/L)*xfac;

arg = 2*K*k*x;
[~,cn,~] = ellipj(arg, m);
u = umin + a*cn.^2; %KdV Solution
T = chi/cg;

% Modify to satisfy zero-mass of Ostrovsky eqn.
integrand = [u; u(1)]; %Compute mean
ubar = 1/(x_max-x_min)*trapz([x;x_max], integrand);

uprime = u - ubar; %zero-mass KdV solution and period
cprime = cg - ubar; 
Tprime = abs(chi/cprime);
%% 
%Only useful for moderate rotation, have to guess T for high rotation
Tstar = EstimateTOstrovskyperiodic(uprime, 15, cprime, M, J, x);
%% Integrations to Visualise Evolution
% Change the input to observe behaviour of minimised values etc.
% Used also to compute momentum if needed as a constraint.

M = 5000; J = 50; %Appear reasonable unless integrating over a long time
[u_lj, u_lj_t, G, dGdT, p_comp, tdata] = ...
    FIOstrovskyEvenxperiodic(u, 1.2*T, cg, M, J, x);
[dGdu, sdata, udata, wdata] = AIOstrovskyEvenxperiodic(u_lj, u_lj_t, Tprime, M, J, x);
%% Waterfall Plot of Evolution
figure;
waterfall(x,tdata,u_lj');
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
%% KdV Initial Values Without Zero Mass Plot
figure;
p1 = plot(x,u_lj(:,1)','color',[0 0.4470 0.7410]);hold on;
plot(x,u_lj(:,2:24)','color',[0 0.4470 0.7410]);hold on;
p2 = plot(x,u_lj(:,25)','color',[0.8500 0.3250 0.0980]);hold on;
plot(x,u_lj(:,26:50)','color',[0.8500 0.3250 0.0980]);
xlabel('x'); ylabel('u(x,t)'); legend([p1,p2],{'descending','ascending'});
%% Comparison of u0 at 0 and End of Time Period Integrated Over
figure
plot(x,u_lj(:,1),x,u_lj(:,J+1));
grid on
xlabel('x'); ylabel('u_0 and u(x,2T)'); legend('u_0^{(min)}','u^{(min)}(x,2T)');
%%
%xmax = @(T) wrapperxmax(uprime, T, cprime, M, J, x);
%Ta = 3.8; Tb = 3.6;
%Tstar = EstimateTBisection(xmax, Ta, Tb);
%%
q = [uprime(1:L/2+1); Tprime]; %Even symmetry
%q = [uprime;Tprime]; %Uneven
[u0min, Tmin, Gmin, niter, info] = ...
    BFGS_OstrovskyFperiodic(q, xfac, M,J,1e-12,30,'Symm',true,'p',0,'A',0);

%% Return Important Values
niter 
info.Gs(1)
Gmin
abs(Tprime+0.0125-Tmin)
1/L*(sum(([u0min;flipud(u0min(2:L/2))]-uprime).^2)) %Even
%1/L*(sum((u0min-uprime).^2))
%% Plot Convergence
figure;
semilogy(0:niter, info.Gs);
grid on;
xlabel('k');
ylabel('G_{tot}(u_0,T)');

%% Return to Fourier Points From Coefficients
% Use this if BFGS_OstrovskyCperiodic was used to return v0min
v0 = fft([u0min;flipud(u0min(2:L/2))]);
%v0 = fft(uprime);
v0(L/2+1) = 0;
v0(8:L-8) = 0;
u01 = real(ifft(v0));
%q = [real(v0(1:L/4)); Tprime];

%%
%function x_max = wrapperxmax(u0, T, c, M, J, x)
%[u_lj, ~, ~, ~, ~, ~, ~] = FIOstrovskyEven(u0, 2*T, c, M, J, x);
%[~, x_i] = max((u_lj(:,J))); x_max = x(x_i);
%end