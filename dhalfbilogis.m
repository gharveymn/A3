function xout=dhalfbilogis(order,xin,c_3,lo,hi,newlo,newhi)
	%Derivative of logit, uc_1 to order 2
	%
	%Author: Gene Harvey
	
		switch nargin
			case 3
				newlo = min(xin);
				newhi = max(xin);
				lo = newlo;
				hi = newhi;
			case 5
				lo = newlo;
				hi = newhi;
			case 6
				error('Too many or too few arguments. Need either 5 or 7, you entered 6.')
		end
	
	c_1 = (exp(c_3/2) + 1)/(exp(c_3/2) - 1);
	
	switch order
		case 1
			xout = (c_1.*c_3.*exp((c_3.*(lo - xin))./(2.*hi - 2.*lo)).*(2.*newhi - 2.*newlo))...
					./((exp((c_3.*(lo - xin))./(2.*hi - 2.*lo)) + 1).^2.*(2.*hi - 2.*lo));
		case 2
			xout = (2.*c_1.*c_3.^2.*exp((2.*c_3.*(lo - xin))./(2.*hi - 2.*lo)).*(2.*newhi - 2.*newlo))...
					./((exp((c_3.*(lo - xin))./(2.*hi - 2.*lo)) + 1).^3.*(2.*hi - 2.*lo).^2) ...
					- (c_1.*c_3.^2.*exp((c_3.*(lo - xin))./(2.*hi - 2.*lo)).*(2.*newhi - 2.*newlo))...
					./((exp((c_3.*(lo - xin))./(2.*hi - 2.*lo)) + 1).^2.*(2.*hi - 2.*lo).^2);
			

	end
	
end