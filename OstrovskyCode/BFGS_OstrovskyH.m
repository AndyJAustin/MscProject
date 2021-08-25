function [u0min, cmin, Tmin, Gmin, niter, info] = BFGS_OstrovskyH(q, b, varargin)
% Implements the BFGS quasi-Newton algorithm to a functional G(u0, c, T)
% that is zero only when the solution u(x,t) of the Ostrovsky equation, with
% initial condition u0(x) = u(x,0), is time periodic with half-period T and
% group speed c. This routine implements the procedure with initial
% condition data (u0) given at Hermite nodes and is therefore most
% appropriate for wavepacket solutions. Peicewise cubic Hermite
% interpolation is responsible for converting the initial condition from
% Hermite points to evenly spaced nodes for pseudo-spectral integration.
% The time taken to perform each iteration is returned to the command
% window.
%
% Input
% -----
% * q: Initial descent variables [u0; c; T] of length d (real). For N Hermite
%      nodes x_n, n=[1,...,N]; u0 consists of first (N-1)/2+1 or all N 
%      values for even symmetry or uneven wavepackets respectively. N must be
%      an odd integer. 
% * b: Scaling parameter x_n = b xtilde_n that allows integration to happen
%      over a narrower or wider range of x. 
% * M: Optional number of time steps. Positive integer. Default value
%      2000. Must be large enough for stability. Computational cost linear
%      in M.
% * J: Optional number of slices of solution matrix u_lm that will form the
%      lower memory solution matrix u_lj j=[1,...,J+1]. I.e. the number of
%      'mileposts' minus one. M mod J = 0, M =/= J. Default value 100.
% * L: Number of Fourier points interpolated into for integration. Must be
%      power of 2 for efficiency. Determines the number of modes for
%      pseudo-spectral integration. Higher L more accurate but more costly
%      than increasing M so ideally should not exceed 1024. Default value
%      512.
% * epsilon: Optional convergence tolerance on the norm(delGk). Positive. 
%      Deafult value 1e-9.
% * MaxIter: Optional maximum number of iterations of BFGS outer loop.
%      Positive integer. Default value 20.
% * Symm: Parameter for symmetry of u0. Logical value true or false. (true
%      means symmetric) Determines symmetry of initial wavepacket data and
%      therefore the representation of q and formulation of objective G and
%      its derivatives.
% * p: Parameter for momentum constraint. Positive or zero; if zero there
%      is no constraint on momentum, if positive, the optimisation
%      constrains the momentum of the solution approximation to p.
% * A: Parameter for constraint on maximum at origin. Positive or zero; if
%      zero there is no constraint on momentum, if positive, the
%      optimisation constrains the momentum of the solution approximation
%      to p.
%
% Optional inputs must keep their order, i.e. if one wanted to change L; M
% and J would have to be entered manually too. Parameters must be preceded
% by their name. I.e. to set even symmetry of u0, one would enter " 'Symm',
% true " direcly after the optional inputs and before parameters p and A
% which in turn must be of the same format, all seperated by a comma.
%
% Output
% ------
%
% * u0min: Optimal value of initial waveform in the appropriate
%       representation. I.e. qmin(1:d-2).
% * cmin: Optimal group speed.
% * Tmin: Optimal half-period.
% * Gmin: Value of objective function at minimum.
% * niter: Number of iterations performed. 
% * info: Struct with fields
%           - Gs: Values of Gk at each iteration.
%           - cs: -||- group speed -||-
%           - Ts: -||- half-period -||-
% 
% One could add more quantities to info such as the wavefunctions at each
% iterate u0_k. We exclude this to limit momory usage and because it can be
% implemented with ease if neccessary.

% Parse Input (Handle any input errors)
% -----------

P = inputParser;
validScalarPos = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosOrZero = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(P,'q');
addRequired(P,'b');
addOptional(P,'M',2000,validScalarPos);                 %Set default values
addOptional(P,'J',100,validScalarPos);
addOptional(P,'L',512,validScalarPos);
addOptional(P,'epsilon',1e-9,validScalarPos);
addOptional(P,'MaxIter',20,validScalarPos);
addParameter(P,'Symm',true,...
              @(x) islogical(x));
addParameter(P,'p',0,validScalarPosOrZero);
addParameter(P,'A',0,validScalarPosOrZero);

parse(P, q, b, varargin{:});

M = P.Results.M;
J = P.Results.J;
L = P.Results.L;
Symm = P.Results.Symm;
pandA = [P.Results.p P.Results.A];
epsilon = P.Results.epsilon;
MaxIter = P.Results.MaxIter;

if P.Results.A == 0 && P.Results.p == 0
    Constraints = 'None';
elseif P.Results.A == 0
    Constraints = 'Momentum';
elseif P.Results.p == 0
    Constraints = 'Amplitude';
else
    Constraints = 'Both';
end
try
    assert(mod(log2(L),1) == 0);
catch
    warning('L must be an integer power of 2 for FFT efficiency; assigning L=512');
end
try 
    assert(mod(M,J) == 0);
catch
    warning('M must be an integer product of J; assigning default values M=2000,J=100');
    M = 2000; J = 100;
end
try 
    assert(M/J >= 1);
catch
    warning('M must be greater than or equal to J; assigning default values M=2000,J=100');
    M = 2000; J = 100;
end

% Set up Functions (For linesearch, as this may only need to compute Gk and
% ----------------  not delGk --> create wrapper functions to calculate Gk 
%                   and delGk individually.)

Phi.G = @(qk) wrapperG(qk);
function Gk = wrapperG(qk) %Nested functions can use global variables
        [Gk, ~] = IntegrateOstrovskyH(qk, b, M, J, L, Symm, Constraints, pandA);     
end
Phi.delG = @(qk) wrapperdelG(qk);
function delGk = wrapperdelG(qk)
        [~, delGk] = IntegrateOstrovskyH(qk, b, M, J, L, Symm, Constraints, pandA);
end

opts.c1 = 1e-4; 
opts.c2 = 0.9;
alpha_max = 1;
linesearch = @(qk, zk) LineSearch(Phi, qk, zk, alpha_max, opts);

% Initialise Descent Algorithm
% ----------------------------
qk = q;
d = length(q);
I = eye(d);
Hk = I;
[Gk, delGk] = IntegrateOstrovskyH(qk, b, M, J, L, Symm, Constraints, pandA);
info.Gs = Gk;
Gkmin1 = Gk+1;
%info.delGs = delGk;
info.cs = qk(d-1);
info.Ts = qk(d);
k = 0;

% Perform BFGS
% ------------
while norm(delGk) > epsilon && k < MaxIter && Gkmin1 - Gk > 0
    tic                                                %time each iteration
    
    % Compute Search Direction
    zk = - Hk*delGk;
    
    % Compute q_kplus1 with Linesearch
    alpha = linesearch(qk, zk);
    qkmin1 = qk;
    qk = qk + alpha*zk;
    
    % Compute Objective, Gradient, sk and yk
    delGkmin1 = delGk;
    Gkmin1 = Gk;
    [Gk, delGk] = IntegrateOstrovskyH(qk, b, M, J, L, Symm, Constraints, pandA);
    
    yk = delGk - delGkmin1;
    sk = qk - qkmin1;
    
    % Scaling Heuristic 
    if k < 1
        Hk = ((yk'*sk)/(yk'*yk))*Hk;
    end
    
    % Compute H_kplus1
    rhok = 1/(yk'*sk);
    Hk = (I - rhok*(sk*yk'))*Hk*(I - rhok*(yk*sk')) + rhok*(sk*sk');
    
    % Store Information
    info.Gs = [info.Gs Gk];
    info.cs = [info.cs qk(d-1)];
    info.Ts = [info.Ts qk(d)];
    
    k = k + 1;
    toc
    
end

% Return Optimal Values
% ---------------------
u0min = qk(1:d-2);
cmin = qk(d-1);
Tmin = qk(d);

Gmin = Gk;
niter = k;

end