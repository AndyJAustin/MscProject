%%           KORTEWEG DE VRIES
%           -------------------

%% Soliton Solutions 
% Only used for understanding and we do not optimise these and instead use
% an x-periodic cnoidal initial conditions.

L = 512;
n = 60;
N = 2*n + 1;
%x_j = herroots(N);
b = 0.3;
x_scaled = x_j./b;
x_max = -1.065*x_scaled(1);
x_min = -x_max; 
xfac = (x_max - x_min)/(2*pi);
xi = (2*pi/L)*(-L/2 : L/2-1)';
x = 0.5*(x_max + x_min) + xfac*xi;

cg = 3;
a0 = 0.1;
k0 = 1.2;
om0 = 1/k0-k0^3;
%cg = -3*k0^2-1/k0^2;
cp = -k0^2 + 1/k0^2;
cgk = -6*k0+2/k0^3;
om2 = 2*k0^3/(12*k0^4+3);
om = om0 + om2*a0^2;
K = a0*sqrt(-om2/cgk);
chi = 2*k0^2/(12*k0^4+3);

D = sqrt(12/a0); %should be 12/a0
u0 = a0*sech(x/D).^2;
T = 30;

%% X-PERIODIC KDV

L = 256;
r1 = 1.6;
r2 = 2;
r3 = 2.1;
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
u0 = umin + a*cn.^2;
T = chi/cg;

%% Integrate Forwards in Time
%Find momentum of exact solution:
M=2000; J=100; %Reasonable numbers
[u_lj, u_lj_t, G, ~, p_comp, tdata] = FIKdVEven(u0prime, Tprime, M, J, x);
%% Initial Values
%Add any slight modification to test effectiveness wrt initial conditions.

u0prime = u0+1e-6*cos(arg); %+1e-3*rand(L,1).*exp(-(.5.*x).^2);%+1e-6*sin(arg);
Tprime = T;
%% Perform Optimisation
% BFGS_KdVC not included as is not effectiev but can be implemented in a
% similar way.

q = [u0prime(1:L/2+1);Tprime]; %Even symmetry
[u0min, Tmin, Gmin, niter, info] = ...
    BFGS_KdVF(q, xfac, M,J,1e-12,30,'Symm',true,'p',p_comp,'A',0);
%%
q = [u0prime; Tprime]; %Uneven symmetry
[u0min, Tmin, Gmin, niter, info] = ...
    BFGS_KdVF(q, xfac, M,J,1e-12,30,'Symm',false,'p',0,'A',0);
%% Return Important Quantities
niter
info.Gs(1)
Gmin
abs(T-Tmin)
%% Plot Convergence
figure;
semilogy(0:niter, info.Gs);
grid on;
xlabel('k');
ylabel('G_{tot}(u_0,T)');
%% Comparison of Initial and Converged u0
figure; %Even
plot(x,[u0min;flipud(u0min(2:L/2))],x,u0prime,'.');
grid on;
xlabel('x');    legend('Converged','Initial');
ylabel('Initial and Converged u_0(x)');
%%
figure; %Uneven
plot(x,u0min,x,u0prime,'.');
grid on;
xlabel('x');    legend('Converged','Initial');
ylabel('Initial and Converged u_0(x)');
%% Used to Smooth Out the Wavefunction
v0 = fft(u0prime);
v0(L/2+1) = 0;
v0(L/4:3*L/4) = 0;
u01 = real(ifft(v0));
q = [real(v0(1:L/4)); Tprime];