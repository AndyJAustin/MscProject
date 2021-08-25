%% Stability of Differentiation Matrices
% herdif.m taken from taken from Weideman and Reddy's MATLAB differentiation
% matrix suite. ACM TOMS, Vol. 26, pp. 465--519 (2000).
% https://appliedmaths.sun.ac.za/~weideman/research/differ.html.

Nmax = 80;
b = 0.1; 
%c = -5;
eig1 = zeros(Nmax,1);
eig2 = zeros(Nmax,1);
eig3 = zeros(Nmax,1);
%eig_plus = zeros(Nmax,1);
for N=2:Nmax
    [~, DM] = herdif(N, 3, b);
    D = DM(:,:,1);
    D2 = DM(:,:,2);
    D3 = DM(:,:,3);
    %D_plus = c*D - D3;
    e_D = eig(D);
    e_D2 = eig(D2);
    e_D3 = eig(D3);
    eig1(N) = max(abs(e_D));
    eig2(N) = max(abs(e_D2));
    eig3(N) = max(abs(e_D3));
    %eig_plus(N) = max(abs(eig(D_plus)));
end

%%
N = 1:Nmax;
figure;
plot(N,eig1,'.', N, b*1.5.*N.^.5);
xlabel('N'); ylabel('\rho'); grid on;
legend('Spectral radius of D_1', '3/2bN^{1/2}');
figure;
plot(N,eig2,'.',N,(b)^2*2.*N);
xlabel('N'); ylabel('\rho'); grid on;
legend('Spectral radius of D_2', '2b^2N');
figure;
plot(N,eig3,'.', N, (b)^3*2.5.*N.^(1.5));
xlabel('N'); ylabel('\rho'); grid on;
legend('Spectral radius of D_3', '5/2b^3N^{3/2}');
