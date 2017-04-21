function xout=halfbilogis(xin,c_3,lo,hi,newlo,newhi)
	%halfbilogis
	%Makes a bijective logistic transform from 0 to 1. That is, f:[0,1]->[0,1] invertible s.t.
	%								/                  c_1              1 \
	%	f(x) = newlo + (2 newhi - 2 newlo) | c_4 + ------------------------- - - |
	%								|          / c_3 (lo - xin) \       2 |
	%								|       exp| -------------- | + 1     |
	%								\          \   2 hi - 2 lo  /         /
	%
	%Author: Gene Harvey.
	
    switch nargin
	    case 2
			lo = min(xin);
			hi = max(xin);
			newlo = lo;
			newhi = hi;
	    case 4
			newlo = lo;
			newhi = hi;
	    case 5
			error('Too many or too few arguments. Need either 2, 4 or 6, you entered 5.')
    end
    
	c_1 = (exp(c_3/2) + 1)/(exp(c_3/2) - 1);
	c_4 = -1/(exp(c_3/2) - 1);
	%c_4 = 1/2 - c_1;
	
	xout = 2.*(newhi-newlo).*(c_1./(1+exp(-c_3*((xin-lo)/(2*(hi-lo))))) + c_4 - 1/2) + newlo;
	
end

