function [u_nj, u_nj_t, G, dGdT, dGdc, momentum, tdata] = FIOstrovskyDiffMat(u_0, T, cg, b, L, M, J)
%For the differentiation matrices approach, we do not need to leave the
%Hermite points except for calculating momentum 
%
% We do not proceed with differentiation matrix approach, and instead opt
% for an integrating factor method, as the computational cost outweigh the
% marginal benefits offered. 
%
% herdif.m taken from taken from Weideman and Reddy's MATLAB differentiation
% matrix suite. ACM TOMS, Vol. 26, pp. 465--519 (2000).
% https://appliedmaths.sun.ac.za/~weideman/research/differ.html.

dt = 2*T/M;
Q = M/J; % number of nodes per slice
N = length(u_0);
u = u_0;

[u0, x] = SplineInterpolant3(u, b, L);
p_integrand = [u0; u0(1)].^2;
momentum = .5*trapz([x;-x(1)], p_integrand);

j = 1;
u_nj = zeros(N, J + 1);
u_nj(:,j) = u;

%Differentiation matrices
[~, DM] = herdif(N, 3, b);
D = DM(:,:,1);
D3 = DM(:,:,3);
D_Nmin1 = D(1:N-1,1:N-1);
dD_Nmin1 = decomposition(D_Nmin1);
u_Nmin1 = u(1:N-1);

% Solve PDE over T
tdata = zeros(J+1,1);
u_nj_t = zeros(N,J+1);
tdata(1) = 0;
for m = 1:M
   t = m*dt;
   %calculate v
   v_Nmin1 = dD_Nmin1\u_Nmin1;
   v = [v_Nmin1;0];
   k1 = (cg*eye(N)*D - diag(u)*D - D3)*u + v;
   
   %calculate v2 (with u+adt/2)
   u_Nmin1 = u_Nmin1 + .5*k1(1:N-1)*dt;
   v_Nmin1 = dD_Nmin1\u_Nmin1;
   v = [v_Nmin1;0];
   k2 = (cg*eye(N)*D - diag(u)*D - D3)*(u+.5*dt*k1) + v;
   
   %calculate v3 with (u+(u+adt/2)dt/2)
   u_Nmin1 = u_Nmin1 + .5*k2(1:N-1)*dt;
   v_Nmin1 = dD_Nmin1\u_Nmin1;
   v = [v_Nmin1;0];
   k3 = (cg*eye(N)*D - diag(u)*D - D3)*(u+.5*dt*k2) + v;
   
   %calculate v4 with (u+(u+(u+adt/2)dt/2)dt/2)
   u_Nmin1 = u_Nmin1 + .5*k3(1:N-1)*dt;
   v_Nmin1 = dD_Nmin1\u_Nmin1;
   v = [v_Nmin1;0];
   k4 = (cg*eye(N)*D - diag(u)*D - D3)*(u+.5*dt*k3) + v;
   
   u = u + dt*(k1 + 2*(k2+k3) + k4)/6;
   
   
   if mod(m, round(Q)) == 0 && m ~= 1
        % add values to u_lj
        j = j + 1;
        u_nj(:, j) = u;

        u_nj_t(:,j) = k1;
        tdata(j) = t; %for visualisation
        if m == M
            u_n2T_x = D*u;
            u_n2T_t = k1;
            u_n2T = u;
        end
        
   elseif m == 1
        u_nj_t(:,j) = k1;
   end
   
   
    
end

%Alternative formulation of objective
[u_l2T, x] = SplineInterpolant3(u_n2T, b, L);
[u_l2T_t, ~] = SplineInterpolant3(u_n2T_t, b, L);
% calculate G and its derivatives wrt T and c
% -- G --
%u_lT_negx = flipud(u_lT);
%G_integrand = (u_lT - u_lT_negx).^2;
%G = .5*trapz(x, G_integrand);
G_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).^2;
G = .5*trapz([x;-x(1)], G_integrand);


% -- dG/dT --
dGdT_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).*[u_l2T_t; u_l2T_t(1)];
dGdT = 2*trapz([x;-x(1)], dGdT_integrand);

% -- dG/dc --
u_l2T_x = SplineInterpolant3(u_n2T_x, b, L);

dGdc_integrand = ([u_l2T; u_l2T(1)] - [u0; u0(1)]).*[u_l2T_x; u_l2T_x(1)];
dGdc = - 2*T*trapz([x;-x(1)], dGdc_integrand);


end
