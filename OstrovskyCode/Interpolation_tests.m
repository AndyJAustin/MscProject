a0 = 1;         %0.3; 1
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

%b = sqrt(2*N)/x_max; %This seems to work reasonably (From Tao Tang SIAM J. ScI. COMPUT.
%Vol. 14, No. 3, pp. 594-606, May 1993)
%b =0.3;
b = 0.38;
x_scaled = x_j./b;
xf_scaled = x_f./b;
x_max = -1.05*x_scaled(1);  % Could try -2L --> +2L
%x_max = -x_scaled(1)-2;
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

a = 2*a0 + 2*chi*a0^2;

alpha = @(x) exp(-.5*x.^2);
alpx = alpha(x.*b);
alpxk = alpha(x_j);
p_Nmin1 = HermiteInterpolation(x_scaled, fk, x, alpxk, alpx);
%p = polint(x_j, fk, x, alpxk, alpx);
%p = polint(x_scaled, fk, x, alpxk, alpx);
p = interp1([x_min;x(3*L/16); x_scaled;x(5*L/16); x_max],[0;0;fk;0;0],x,'pchip');
MSE_Hermite = sum((p_Nmin1-u0).^2)/L;
MSE_polint = sum((p - u0).^2)/L;

T = 0.8; 
M = round(T/((0.4/L^2) * xfac^2));
J = 100;
%%
b2 = 0.15;
x_scaled2 = x_j./b2;
alpx2 = alpha(x_scaled2.*b);
p2 = polint(x_scaled, fk, x_scaled2, alpxk, alpx2);


%%
%fk = zeros(N,1);
%fk(n+1) = a0;
%fk(n+2:N) = u0;
%fk(1:n) = flipud(u0);
%alpha = @(z) exp(-z.^2/2);
%alpxk = alpha(x_j);
%alpx = alpha(x.*b);
%u_0 = polint(x_scaled, fk, x, alpxk, alpx);

%Interpolate back to Hermite points
%x_j = herroots(N);
%w = real(ifft(w_hat));
%w = w + flipud(w);
%w_n = fourint(w, x_j/(be*xfac));
%Return dG/du_0
%n = (N-1)/2;
%dGdu = w_n(n+2:N);

%%
tic
hermiteH(81,x);
toc
tic
polint(x_scaled, fk, x, alpxk, alpx);
toc
