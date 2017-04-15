function xout=dbilogit(order,xin,c_3,lo,hi,newlo,newhi)
	%Derivative of logit, up to order 2
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
	c_4 = -1/(exp(c_3/2) - 1);
	
	switch order
		case 1
			xout = -(c_1.*(newhi - newlo))./(c_3.*(c_1./(c_4 + (lo - xin)./(hi - lo)) + 1)...
				.*(hi - lo).*(c_4 + (lo - xin)./(hi - lo)).^2);
		case 2
			xout = (c_1.^2*(newhi - newlo))./(c_3.*(c_1./(c_4 + (lo - xin)./(hi - lo)) + 1).^2 ...
					.*(hi - lo).^2.*(c_4 + (lo - xin)./(hi - lo)).^4) ...
					- (2.*c_1.*(newhi - newlo))./(c_3.*(c_1./(c_4 + (lo - xin)./(hi - lo)) + 1)...
					.*(hi - lo).^2.*(c_4 + (lo - xin)./(hi - lo)).^3);
end

