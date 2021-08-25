%%          OSTROVSKY
%           ----------
%Second order packet defined by a0 and k0 (must guess for T)
a0 = 1;         %0.3;
k0 = 0.85;    %0.93; 0.81
om0 = 1/k0-k0^3;            cg = -3*k0^2-1/k0^2;
cp = -k0^2 + 1/k0^2;        cgk = -6*k0+2/k0^3;
om2 = 2*k0^3/(12*k0^4+3);   om = om0 + om2*a0^2;
K = a0*sqrt(-om2/cgk);      chi = 2*k0^2/(12*k0^4+3);

L = 512;
n = 60;
N = 2*n + 1;
x_j = herroots(N);
x_f = x_j(1:n+1);

%b =0.3;
b = 0.38;
x_scaled = x_j./b;
xf_scaled = x_f./b;
x_max = -1.05*x_scaled(1);
x_min = -x_max; 
xfac = (x_max - x_min)/(2*pi);
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = 0.5*(x_max + x_min) + xfac*xi;

% Second order packet
u0 = 2.*a0.*sech(K*x).*cos(k0*x)+2.*(a0^2).*chi.*((sech(K*x)).^2).*cos(2*k0*x);
fk = 2.*a0.*sech(K*x_scaled).*cos(k0*x_scaled)+2.*(a0^2).*chi.*((sech(K*x_scaled)).^2).*cos(2*k0*x_scaled);
f = 2.*a0.*sech(K*xf_scaled).*cos(k0*xf_scaled)+2.*(a0^2).*chi.*((sech(K*xf_scaled)).^2).*cos(2*k0*xf_scaled);
fk2 = 2.*a0.*sech(K*x_j).*cos(k0*x_j)+2.*(a0^2).*chi.*((sech(K*x_j)).^2).*cos(2*k0*x_j);

v0 = fft(u0);
v0(L/2+1) = 0; %Nyquist frequency
v0(L/8:7*L/8) = 0;
u01 = real(ifft(v0));

T = 0.8; 
%M = round(T/((0.4/L^2) * xfac^2));
J = 100;
%%
Tstar = EstimateTOstrovsky(u0, T, 2000, J, x);
%% Call to Forward and Adjoint Integration Functions
% Change the input values as desired, u_lj, tdata can be used to plot
% evolution. momentum p_comp computed this way.

M = 2000; J = 100;
[u_lj, u_lj_t, G, dGdT, dGdc, p_comp, tdata] = ...
    FIOstrovskyEven(u0, 2*T, cg, M, J, x);
[dGdu, sdata, udata, wdata] = AIOstrovskyEven(u_lj, u_lj_t, Tprime, cg, M, J, x);
%% Waterfall plot
figure;
waterfall(x,tdata,u_lj');
xlabel('x'); ylabel('t'); zlabel('u(x,t)');
%% u0 Evolution frames
figure
plot(x,u_lj(:,1),x,u_lj(:,J+1));
grid on
xlabel('x'); ylabel('u_0 and u(x,2T)'); legend('u_0^{(min)}','u^{(min)}(x,2T)');
%% u0 Comparison
figure
%plot(x,[u0min;flipud(u0min(2:L/2))],x,u0,'-.'); %If even
plot(x,u0min,x,u0,'-.');
grid on
xlabel('x'); ylabel('u_0^{(0)} and u_0^{(min)}'); legend('u_0^{(min)}','u_0^{(0)}');
%% Optimisation w. Fourier Nodes
% Change optimisation parameters if neccessary
% Set up q0 first
q = [u0(1:L/2+1); cg; 0.8238]; %Even
%q = [u0;cg;Tprime]; %Uneven
[u0min,cmin, Tmin, Gmin, niter, info] = ...
    BFGS_OstrovskyF(q, xfac, 2000,J,1e-12,50,'Symm',true,'p',0,'A',0);

%% Optimisation w. Hermite Nodes
q = [f; cg; 0.8326];
[u0_min,c_min, T_min, Gmin, niter, info] = ...
    BFGS_OstrovskyH(q, xfac,2000,100,L,1e-12,100,'Symm',true,'p',0,'A',0);

%% Convert Back from Hermite Nodes
[u0min, x, cmin, Tmin] = H2F([u0_min;c_min;T_min], b, L, true);
%% Return Important Quantities
niter %Return values
info.Gs(1)
Gmin
abs(Tstar-Tmin)
1/L*(sum((u0min-u0).^2)) %Even
%1/L*(sum((u0min-uprime).^2))
%% Info Plots
figure;
semilogy(0:niter, info.Gs);
grid on;
xlabel('k');
ylabel('G_{tot}(u_0,T)');
figure;
plot(0:niter, info.Ts);
grid on;
xlabel('k');
ylabel('T_k');
figure;
plot(0:niter, info.cs);
grid on;
xlabel('k');
ylabel('c_k');