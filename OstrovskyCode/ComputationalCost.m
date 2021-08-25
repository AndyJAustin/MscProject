%% Computational Cost of Integration
a0 = 1;         %0.3; 1
k0 = 0.85;    %0.93; 0.81
om0 = 1/k0-k0^3;
cg = -3*k0^2-1/k0^2;
cp = -k0^2 + 1/k0^2;
cgk = -6*k0+2/k0^3;
om2 = 2*k0^3/(12*k0^4+3);
chi = 2*k0^2/(12*k0^4+3);
a = 2*a0 + 2*chi*a0^2;
om = om0 + om2*a^2;
K = a0*sqrt(-om2/cgk);

%T = 0.125*pi/om;
T = 0.848;
M = 2000;
J = 100;
b = 0.3;
%% Compare Hermite Differentiation Matrices with Integrating Factor Method.
% Outputs relevant data to the command window.
Ns = [31,41,51,61,71,81,91,101,111,121,131,141,151,161,171,181];
Ls = [32,64,128,256,512,1024];

Nts = zeros(length(Ns),1);
for i = 1:length(Ns)
    N = Ns(i)
    x_j = herroots(N);
    x_scaled = x_j./b;
    f = 2.*a0.*sech(K*x_scaled).*cos(k0*x_scaled)+...
        2.*(a0^2).*chi.*((sech(K*x_scaled)).^2).*cos(2*k0*x_scaled);
    tic
    [u_nj, u_nj_t, G, dGdT, dGdc, momentum, tdata] = ...
        FIOstrovskyDiffMat(f, T, cg, b, 512, M, J); %L only for calculating derivatives
    toc
end
x_max = 50; x_min = -x_max;
xfac = (x_max - x_min)/(2*pi);
Lts = zeros(length(Ls),1);
for i = 1:length(Ls)
    L = Ls(i)
    xi = (2*pi/L)*(-L/2 : L/2-1)';
    x = xfac*xi;
    u0 = 2.*a0.*sech(K*x).*cos(k0*x)+...
        2.*(a0^2).*chi.*((sech(K*x)).^2).*cos(2*k0*x);
    tic
    [u_lj, u_lj_t, G, dGdT, dGdc, momentum, tdata] = ...
        FIOstrovskyOdd(u0, T, cg, M, J, x);
    toc
end
%%
figure;
semilogy(Ns,Nts,'.',Ls,Lts,'.'); grid on;
xlabel('N and L'); ylabel('Computation time (s)');
legend('Hermite','Foruier');