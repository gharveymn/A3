function xout=bilogit(xin,c_3,lo,hi,newlo,newhi)
	%Inverse of bilogis
	%
	%Author: Gene Harvey
	
    switch nargin
	    case 4
		   newlo = lo;
		   newhi = hi;
	    case 5
		   error('Too many or too few arguments. Need either 4 or 6, you entered 5.')
    end
	
	c_1 = (exp(c_3/2) + 1)/(exp(c_3/2) - 1);
	c_4 = -1/(exp(c_3/2) - 1);
	
	xout = newlo - (log(- c_1./(c_4 + (lo - xin)./(hi - lo)) - 1)./c_3 - 1/2).*(newhi - newlo);
	
end

