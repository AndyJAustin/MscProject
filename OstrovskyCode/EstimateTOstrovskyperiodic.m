function Tstar = EstimateTOstrovskyperiodic(u0, T0, M, J, x)

if c<0
    sign = -1;
else
    sign = 1;
end

[u_lj, ~, ~, ~, ~, ~] = FIOstrovskyEvenxperiodic(u0, 2*T0, M, J, x);
[~,xindex] = max(u_lj(:,J));
xmax = x(xindex);

k = 0; 
Tk = T0;
while k < 300 && abs(xmax) > 1e-9
    if xmax < 0
        Tk = Tk + sign*xmax/(1e3*abs(c));
    else
        Tk = Tk - sign*xmax/(1e3*abs(c));
    end
    [u_lj, ~, ~, ~, ~, ~] = FIOstrovskyEvenxperiodic(u0, 2*Tk, M, J, x);
    [~,xindex] = max(u_lj(:,J));
    xmax = x(xindex);
    k = k+1;   
end
Tstar = Tk; k
end