function [Gk, delGk] = IntegrateOstrovskyH(qk, b, M, J, L, Symm, Constraints, pandA)
% Performs forward and backward integrations and returns the value of G and
% its gradient at the current value of qk. Takes care of the
% transformations between Hermite points and Fourier points along with the
% symmetry of u0 and any compination of constraints or lack thereof.
%
% Input
% -----
% * qk: Descent variables [u0; c; T] at current iterate.
% * b: Scaling parameter x = b xtilde.
% * M: Number of time-steps.
% * J: Number of slices to approximate u(x,s) in Adjoint integration.
% * L: Number of Fourier points/modes.
% * Symm: Even symmetry of u0 (true/false).
% * Constraints: Constraints on momentum flux and mx value at origin.
% * pandA: A 2X1 array of the values of momentum flux and amp.
%
% Output
% ------
% * Gk: Value of the objective at qk.
% * delGk: Value of the gradient at qk.

% Convert to Fourier Points
% -------------------------
[u0, x, c, T] = H2F(qk, b, L, Symm);

% Perform Forward and Adjoint Integrations
% ----------------------------------------
if isequal(Symm, true)
    [u_lj, u_lj_t, G, dGdT, dGdc, p_comp, ~] = FIOstrovskyEven(u0, T, c, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyEven(u_lj, u_lj_t, T, c, M, J, x);
    N = 2*(length(qk) - 3) + 1;
else
    [u_lj, u_lj_t, G, dGdT, dGdc, p_comp, ~] = FIOstrovskyUneven(u0, T, c, M, J, x);
    [dGduF, ~, ~] = AIOstrovskyUneven(u_lj, u_lj_t, T, c, M, J, x);
    N = length(qk) - 2;
end


% Compute G_tot and its Gradient wrt qk            (Depends on constraints)
% -------------------------------------
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

% Convert Back to Hermite Points Representation
% ---------------------------------------------
[dGduH] = F2H(dGduF, b, N, Symm);
delGk = [dGduH; dGdc; dGdT];

end