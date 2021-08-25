function [v0min, cmin, Tmin, Gmin, niter, info] = BFGS_OstrovskyC(q, xfac, varargin)
% Implements the BFGS quasi-Newton algorithm to a functional G(u0, c, T)
% that is zero only when the solution u(x,t) of the Ostrovsky equation, with
% initial condition u0(x) = u(x,0), is time periodic with half-period T and
% group speed c. This routine implements the procedure with initial
% condition data (u0) represented by Fourier coefficients c_k = a_k + b_k i. 
% As u0 is real, the only modification to q is if the wavefucntion is
% constrained to even symmetry about x=0, in this case we ignore the
% complex coefficients b_k which alters d; the size of q.
%
% The time taken to perform each iteration is
% returned to the command window.
%
% Input
% -----
% * q: Initial descent variables [v0; c; T] of length d (real). Where v0 is
%      the Fourier transform of u0 (in k-space) For L Fourier points (a
%      power of 2 for FFT efficiency), d takes the value 
%      of L/2+1 or L/4 + 1 for a wavefunction that is constrained by even
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
% * v0min: Optimal value of initial waveform in Fourier k-space in the
%          appropriate representation. I.e. qmin(1:d-2).
% * cmin: Optimal group speed.
% * Tmin: Optimal half-period.
% * Gmin: Value of objective function at minimum.
% * niter: Number of iterations performed. 
% * info: Struct with fields
%           - Gs: Values of Gk at each iteration.
%           - cs: -||- group speed -||-
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
addParameter(P,'Symm',true,...          %Symmetry (Even/Uneven)
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

Phi.G = @(qk) wrapperG(qk);
function Gk = wrapperG(qk) %Nested functions can use global variables
        [Gk, ~] = IntegrateOstrovskyC(qk, xfac, M, J, Symm, Constraints, pandA);     
end
Phi.delG = @(qk) wrapperdelG(qk);
function delGk = wrapperdelG(qk)
        [~, delGk] = IntegrateOstrovskyC(qk, xfac, M, J, Symm, Constraints, pandA);
end

opts.c1 = 1e-4; %These appear to be optimal for the Fourier coefficient rep.
opts.c2 = 0.45;
alpha_max = 1;
linesearch = @(qk, zk) LineSearch(Phi, qk, zk, alpha_max, opts);

%---------------------BFGS--------------------
%---------------------------------------------
qk = q;
d = length(q);
I = eye(d);
Hk = I;
[Gk, delGk] = IntegrateOstrovskyC(qk, xfac, M, J, Symm, Constraints, pandA);
info.Gs = Gk;
%info.delGs = delGk;
info.cs = qk(d-1);
info.Ts = qk(d);
k = 0;
while norm(delGk) > epsilon && k < MaxIter
    tic
    %compute search direction
    zk = - Hk*delGk;
    
    %compute x_kplus1 with linesearch
    alpha = linesearch(qk, zk);
    qkmin1 = qk;
    qk = qk + alpha*zk;
    
    %compute objective, delG_kplus1, sk and yk
    delGkmin1 = delGk;
    [Gk, delGk] = IntegrateOstrovskyC(qk, xfac, M, J, L, Symm, Constraints, pandA);
    
    yk = delGk - delGkmin1;
    sk = qk - qkmin1;
    
    if k < 1
        Hk = ((yk'*sk)/(yk'*yk))*Hk;
    end
    %compute H_kplus1
    rhok = 1/(yk'*sk);
    Hk = (I - rhok*(sk*yk'))*Hk*(I - rhok*(yk*sk')) + rhok*(sk*sk');
    
    info.Gs = [info.Gs Gk];
    info.cs = [info.cs qk(d-1)];
    info.Ts = [info.Ts qk(d)];
    
    k = k + 1;
    toc
end %(return values)

v0min = qk(1:d-2);
cmin = qk(d-1);
Tmin = qk(d);

Gmin = Gk;
niter = k;

end