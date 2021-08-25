function alphastar = zoom(phi, alphalo, alphahi, c1, c2)
% Refer to Algorithm 3.3 from Nocedal & Wright, Numerical Optimization,
% Springer 2000
% 
% Input
% -----
% * Phi: Struct with fields
%       - G: function handler for objective as a function of alpha at
%       current qk
%       - delG: gradient handler for gradient of objective
% * alphalo: lower boundary of the interval
% * alphahi: upper boundary of the interval
% * c1: constant in sufficient decrease condition 
%           G(qk + alphak zk) > G(qk) + c1 alphak (delGk'*zk)
%           default 1e-4
% * c2: constant in strong curvature condition 
%           |delG(qk + alphak zk)'*zk| <= c2 |df(qk)'*zk| 
%           default 0.9
%         0 < c1 < c2 < 1
%
% Output
% ------
% * alphastar: optimal step length in interval
%
% Generates an iterate alphak between alphahi and alphalo, then replaces on
% of these end points such that the following are still satisfied.
% a) The interval bounded by alphalo and alphahi contains step lengths that
% satisfy the strong Wolfe conditions
% b) alphalo is, among all step lengths generated so far and satisfying the
% sufficient decrease condition
% c) alphahi is chosen so that delG'(alphalo)(alphahi-alphalo)<0

k = 1;
stop = false;
maxIter = 50;
while k < maxIter && stop == false
    % Interpolate using bisection to find trial step length between
    % alphalow and alphahigh
    alphak = .5*(alphalo + alphahi);
      
    phik = phi.G(alphak);
    
    if abs(alphahi - alphalo) < eps
      alphastar = alphak;
      stop = true;
    end
    
    if phik > phi.G(0) + c1*alphak*phi.delG(0) || phi.G(alphak) >= phi.G(alphalo)
      alphahi = alphak;
      
    else 
        delphi_j = phi.delG(alphak);
            
        if abs(delphi_j) <= -c2*phi.delG(0) % sufficient decrease condition satisfied
            alphastar = alphak;
            stop = true;
        elseif (alphahi - alphalo)*delphi_j >= 0
            alphahi = alphalo;
        end
        alphalo = alphak;
    end
    k = k + 1;

end
if k == maxIter
    alphastar = alphak;
end
end
            
