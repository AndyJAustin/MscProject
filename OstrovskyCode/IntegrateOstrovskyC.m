function [Gk, delGk] = IntegrateOstrovskyC(qk, xfac, M, J, Symm, Constraints, pandA)
% Performs forward and backward integrations and returns the value of G and
% its gradient at the current value of qk. Takes care of the symmetry of u0
% and any compination of constraints or lack thereof, and Determines L from qk.
%
% Input
% -----
% * qk: Descent variables [u0; c; T] at current iterate.
% * b: Scaling parameter x = b xtilde.
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
v0 = qk(1:d-2);
c = qk(d-1); T = qk(d);
[u0, x] = C2F(v0, xfac, Symm);

% ----- Perform Integrations to Compute Derivatives -----

if isequal(Symm, true)
    [u_lj, u_lj_t, G, dGdT, dGdc, p_comp, ~] = FIOstrovskyEven(u0, T, c, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyEven(u_lj, u_lj_t, T, c, M, J, x);
else
    [u_lj, u_lj_t, G, dGdT, dGdc, p_comp, ~] = FIOstrovskyUneven(u0, T, c, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyUneven(u_lj, u_lj_t, T, c, M, J, x);
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

[dGduC] = F2C(dGduF, Symm);
d = length(qk);
delGk = zeros(d,1);
if isequal(Symm, true)
    delGk(1:d-2) = dGduC(1:L/4);
    delGk(d-1) = dGdc;
    delGk(d) = dGdT;
else
    delGk(1:d-2) = dGduC;
    delGk(d-1) = dGdc;
    delGk(d) = dGdT;
end

end