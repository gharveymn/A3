function xout=dbilogis(xin,c_3,order,newlo,newhi,lo,hi)
	%Derivative of logis, up to order 2
	%
	%Author: Gene Harvey
	
		switch nargin
			case 3
				lo = min(xin);
				hi = max(xin);
				newlo = lo;
				newhi = hi;
			case 5
				newlo = lo;
				newhi = hi;
			case 6
				error('Too many or too few arguments. Need either 5 or 7, you entered 6.')
		end
	
	c_1 = (exp(c_3/2) + 1)/(exp(c_3/2) - 1);
	
	switch order
		case 1
			xout = (c_1.*c_3.*exp(c_3.*((lo - xin)./(hi - lo) + 1/2)).*(newhi - newlo))...
					./((hi - lo).*(exp(c_3.*((lo - xin)./(hi - lo) + 1/2)) + 1).^2);
		case 2
			xout = (2*c_1*c_3^2*exp(2*c_3*((lo - xin)./(hi - lo) + 1/2)).*(newhi - newlo))...
					./((hi - lo).^2*(exp(c_3*((lo - xin)./(hi - lo) + 1/2)) + 1).^3)...
					- (c_1*c_3^2*exp(c_3*((lo - xin)./(hi - lo) + 1/2)).*(newhi - newlo))...
					./((hi - lo).^2*(exp(c_3*((lo - xin)./(hi - lo) + 1/2)) + 1).^2);
	end
end