function [u0min, Tmin, Gmin, niter, info] = BFGS_KdVF(q, xfac, varargin)
% Implements the BFGS quasi-Newton algorithm to a functional G(u0, T)
% that is zero only when the solution u(x,t) of the KdV equation, with
% initial condition u0(x) = u(x,0), is time periodic with half-period T and
% group speed c. This routine implements the procedure with initial
% condition data (u0) represented by evenly spaced Fourier points. There is
% no need for any conversion to an alternte representation as all of the
% computation occurs within this domain. The modification that us made to
% u0 before integrating is the reflection about x=0 if an even symmetry
% constraint is chosen. The time taken to perform each iteration is
% returned to the command window. As we are in a
% periodic x domain, we perform the integration from a fixed frame and do
% not include a translation term of the KdV eqaution, negating any need
% for a speed variable c. (We are focused on the cnoidal solutions not
% soliton solutions as these create weak minimisers) The time taken to
% perform each iteration is returned to the command window. 
%
% Input
% -----
% * q: Initial descent variables [u0; T] of length d (real). For L
%      Fourier points (a power of 2 for FFT efficiency), d takes the value
%      of L/2+2 or L + 2 for a wavefunction that is constrained by even
%      symmetry about the origin, and not so, respectively. L is computed
%      from the size of q and the decision to implement an even symmetry or
%      not.
% * xfac: Scaling parameter xfac=(xmax-xmin)/2pi that allows integration to happen
%      over a narrower or wider range of x. 
% * M: Optional number of time steps. Positive integer. Default value
%      2000. Must be large enough for stability. Computational cost linear
%      in M.
% * J: Optional number of slices of solution matrix u_lm that will form the
%      lower memory solution matrix u_lj j=[1,...,J+1]. I.e. the number of
%      'mileposts' minus one. M mod J = 0, M =/= J. Default value 100.
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
% * Tmin: Optimal half-period.
% * Gmin: Value of objective function at minimum.
% * niter: Number of iterations performed. 
% * info: Struct with fields
%           - Gs: Values of Gk at each iteration.
%           - Ts: -||- half-period -||-

%Parse input
P = inputParser;
validScalarPos = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalarPosOrZero = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
addRequired(P,'q');
addRequired(P,'xfac');
addOptional(P,'M',2000,validScalarPos); %Setting default values
addOptional(P,'J',100,validScalarPos);
addOptional(P,'epsilon',1e-9,validScalarPos);
addOptional(P,'MaxIter',20,validScalarPos);
addParameter(P,'Symm',true,...
              @(x) islogical(x));
addParameter(P,'p',0,validScalarPosOrZero);
addParameter(P,'A',0,validScalarPosOrZero);

parse(P, q, xfac, varargin{:});

M = P.Results.M;
J = P.Results.J;
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

Phi.G = @(qk) wrapperG(qk);
function Gk = wrapperG(qk) %Nested functions can use global variables
        [Gk, ~] = IntegrateKdVF(qk, xfac, M, J, Symm, Constraints, pandA);
end
Phi.delG = @(qk) wrapperdelG(qk);
function delGk = wrapperdelG(qk)
        [~, delGk] = IntegrateKdVF(qk, xfac, M, J, Symm, Constraints, pandA);
end

opts.c1 = 1e-4;   
opts.c2 = 0.1;    % Normally 0.9 but this is ineffective
alpha_max = 1;
linesearch = @(qk, zk) LineSearchKdV(Phi, qk, zk, alpha_max, opts);

qk = q;
d = length(q);
I = eye(d);
Hk = I;
[Gk, delGk] = IntegrateKdVF(qk, xfac, M, J, Symm, Constraints, pandA);
info.Gs = Gk;
%info.delGs = delGk;
info.Ts = qk(d);
k = 0;
Gkmin1 = Gk+1;
while norm(delGk) > epsilon && k < MaxIter && Gkmin1-Gk > 0
    tic
    %compute search direction
    zk = - Hk*delGk;
    
    %compute x_kplus1 with linesearch
    alpha = linesearch(qk, zk);
    qkmin1 = qk;
    qk = qk + alpha*zk;
    
    %compute objective, delG_kplus1, sk and yk
    Gkmin1 = Gk;
    delGkmin1 = delGk;
    [Gk, delGk] = IntegrateKdVF(qk, xfac, M, J, Symm, Constraints, pandA);
    
    yk = delGk - delGkmin1;
    sk = qk - qkmin1;
    
    if k < 1
        Hk = ((yk'*sk)/(yk'*yk))*Hk;
    end
    %compute H_kplus1
    rhok = 1/(yk'*sk);
    Hk = (I - rhok*(sk*yk'))*Hk*(I - rhok*(yk*sk')) + rhok*(sk*sk');
    
    info.Gs = [info.Gs Gk];
    info.Ts = [info.Ts qk(d)];
    
    k = k + 1;
    toc
end %(return values)

u0min = qk(1:d-1);
Tmin = qk(d);

Gmin = Gk;
niter = k;

end