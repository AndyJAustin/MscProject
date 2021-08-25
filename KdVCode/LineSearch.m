function alphastar = LineSearch(Phi, qk, zk, alphamax, varargin)
% Refer to Algorithm 3.2 from Nocedal & Wright, Numerical Optimization,
% Springer 2000
%
% Input
% -----
% * Phi: Struct with fields
%       - G: function handler for objective as a function of q_k
%       - delG: gradient handler for gradient of objective
% * qk: Current iterate
% * zk: Current descent vector
% * alphamax: Maximum step length 
% * opts: Optional paremeters
%       - c1: constant in sufficient decrease condition 
%           G(qk + alphak zk) > G(qk) + c1 alphak (delGk'*zk)
%           default 1e-4
%       - c2: constant in strong curvature condition 
%           |delG(qk + alphak zk)'*zk| <= c2 |df(qk)'*zk| 
%           default 0.9
%         0 < c1 < c2 < 1
%
% Output
% ------
% alphastar: Optimal step length
%
% Procedure uses the knowledge that the interval (alphai_min1, alpha_i)
% contain step lengths that satisfy the strong Wolfe conditions if one of
% the following conditions holds:
%
% i) alphai violates the sufficient decrease condtition
% ii) G(alphai) >= G(alphai_min1)
% iii) delG(alphai) >= 0


P = inputParser;
validScalarZerotoOne = @(x) isnumeric(x) && isscalar(x) ... 
    && (x > 0) && (x < 1);
% 0 < c1 < c2 < 1
addRequired(P,'Phi');
addRequired(P,'qk');  
addRequired(P,'zk');
addRequired(P,'alphamax');
addOptional(P,'c1',1e-4,validScalarZerotoOne);
addOptional(P,'c2',0.9,validScalarZerotoOne);
%Exception for c2 < c1

parse(P, Phi, qk, zk, alphamax, varargin{:});
c1 = P.Results.c1;
c2 = P.Results.c2;
try
    assert(c1 < c2 && 0 < c1 && c1 < 1 && 0 < c2 && c2 < 1);
catch
    warning('c1 must be less than c1 and both must lie in (0,1)');
end

% Set up functions that only take alpha as an input
phi.G = @(alpha) Phi.G(qk + alpha*zk);
phi.delG = @(alpha) (Phi.delG(qk + alpha*zk)')*zk;

% Initialise parameters and values
factor = 1.5; 
alphai_min1 = 0;
Gi_zero = phi.G(0);
Gi_min1 = Gi_zero;
delGi_zero = phi.delG(0);
alphai = 0.1*alphamax;
alphastar = 0;
i = 2;
maxIter = 10;
stop = false;

while (i < maxIter && stop == false)
    
    Gi = phi.G(alphai);
    delGi = phi.delG(alphai);
    if(Gi > Gi_zero + c1*alphai*delGi_zero || (Gi >= Gi_min1 && i > 2))
        alphastar = zoom(phi, alphai_min1, alphai, c1, c2);
        stop = true;
    elseif(abs(delGi) <= -c2*delGi_zero)
        alphastar = alphai;
        stop = true;
    elseif(delGi >= 0)
        alphastar = zoom(phi, alphai, alphai_min1, c1, c2);
        stop = true;
    end
    alphai_min1 = alphai;
    alphai = min(factor*alphai, alphamax);  %Geometric
    %alphai = (factor*alphai + alpha_max)/(FACT + 1); %Weighted mean
    %alphai = 0.5*(alphai + alpha_max);  %Bisection
    
    Gi_min1 = Gi;
    
    i = i + 1;
end  
