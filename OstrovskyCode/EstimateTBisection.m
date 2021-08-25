function Tstar = EstimateTBisection(xmax, Ta, Tb)

if xmax(Ta)*xmax(Tb)>0 
    return 
else
    k = 0;
    Tmid = (Ta + Tb)/2;
    error = abs(xmax(Tmid));
    while k < 100 && error > 1e-6
        if xmax(Tmid)*xmax(Ta) < 0 
            Tb = Tmid;
        else
            Ta = Tmid;    
        end
        Tmid = (Ta + Tb)/2; 
        error = abs(xmax(Tmid));
        k = k+1;
    end
    Tstar = Tmid;
end