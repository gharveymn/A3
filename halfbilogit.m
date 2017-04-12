function xout=halfbilogit(xin,c_3,lo,hi,newlo,newhi)
	%Inverse of halfbilogis.
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
			error('Too many or too few arguments. Need either 4 or 6, you entered 5.')
    end
    
	c_1 = (exp(c_3/2) + 1)/(exp(c_3/2) - 1);
	c_4 = -1/(exp(c_3/2) - 1);
	%c_4 = 1/2 - c_1;
	
	xout = newlo - (log((2.*c_1.*hi - 2.*c_1.*lo)...
			./(hi - 2.*lo + xin - 2.*c_4.*hi + 2.*c_4.*lo) - 1).*(2.*newhi - 2.*newlo))./c_3;
	
end


