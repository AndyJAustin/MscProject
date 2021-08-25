function [Gk, delGk] = IntegrateOstrovskyFperiodic(qk, xfac, M, J, Symm, Constraints, pandA)
% Performs forward and backward integrations for a PERIODIC DOMAIN and
% returns the value of G and its gradient at the current value of qk. Takes
% care of the symmetry of u0 and any compination of constraints or lack
% thereof, and Determines L from qk. Epsilon (the rotation parameter) can
% be controlled within the FI and AI integration sub-rotines.
%
% Input
% -----
% * qk: Descent variables [u0; T] at current iterate.
% * xfac: Scaling parameter xfac = (xmax-xmin)/2pi.
% * M: Number of time-steps.
% * J: Number of slices to approximate u(x,s) in Adjoint integration.
% * Symm: Even symmetry of u0 (true/false).
% * Constraints: Constraints on momentum flux and mx value at origin.
% * pandA: A 2X1 array of the values of momentum flux and amp.
%
% Output
% ------
% * Gk: Value of the objective at qk.
% * delGk: Value of the gradient at qk.

d = length(qk);
T = qk(d);

% ----- Perform Integrations to Compute Derivatives -----

if isequal(Symm, true)
    L = 2*(d-2);
    u0 = zeros(L,1);
    xi = (2*pi/L)*(-L/2 : L/2-1)';
    x = xfac*xi;
    u0(1:L/2+1) = qk(1:d-1);
    u0(L/2+2:L) = flipud(qk(2:d-2));
    [u_lj, u_lj_t, G, dGdT, p_comp, ~] = FIOstrovskyEvenxperiodic(u0, T, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyEvenxperiodic(u_lj, u_lj_t, T, M, J, x);
else
    L = d-1;
    u0 = zeros(L,1);
    xi = (2*pi/L)*(-L/2 : L/2-1)';
    x = xfac*xi;
    u0(1:L) = qk(1:d-1);
    [u_lj, u_lj_t, G, dGdT, p_comp, ~] = FIOstrovskyUnevenxperiodic(u0, T, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyUnevenxperiodic(u0, u_lj, u_lj_t, T, M, J, x);
end

% ----- Account for any Constraints -----

switch Constraints
    case 'None'
        Gk = G;
    case 'Momentum'
        Gp = .5*(p_comp - pandA(1))^2;
        Gk = G + Gp;
        dGpdu = (p_comp - pandA(1))*u0;
        dGduF = dGduF + dGpdu;
    case 'Amplitude'
        j = length(u0)/2 + 1;
        Gp = .5*(u0(j) - pandA(2))^2;
        Gk = G + Gp;
        dGduF(j) = u0(j) - pandA(2);
    case 'Both'
        j = length(u0)/2 + 1;
        Gp = .5*(p_comp - pandA(1))^2 + .5*(u0(j) - pandA(2))^2;
        Gk = G + Gp;
        dGpdu = (p_comp - pandA(1))*u0;
        dGduF(j) = u0(j) - pandA(2);
        dGduF = dGduF + dGpdu;
end

% ----- Convert back to q representation -----

delGk = zeros(d,1);
if isequal(Symm, true)
    delGk(1:d-1) = dGduF(1:L/2+1); delGk(d) = dGdT;
else
    delGk(1:d-1) = dGduF; delGk(d) = dGdT;
end

end