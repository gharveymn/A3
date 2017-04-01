function [omg,eigVects,r]=getDataFull(numKVals,kVals,h,interval,numEigs)
		
	%Build Matrices
	r = (interval(1):h:interval(2))';
	sz = size(r,1);
	rInv = 1./r;
	rInvM = diag(rInv);
	ident = eye(sz);
	
	deriv = toeplitz([0;-1;zeros(sz-2,1)],[0,1,zeros(1,sz-2)]);
	deriv = deriv./(2*h);
	deriv(1,1) = -1/h;
	deriv(1,2) = 1/h;
	deriv(sz,sz) = 1/h;
	deriv(sz,sz-1) = -1/h;
	deriv2 = deriv^2;
	deriv2(2,:) = deriv2(3,:);
	deriv2(1,:) = deriv2(2,:);
	deriv2(end-1,:) = deriv2(end-2,:);
	deriv2(end,:) = deriv2(end-1,:);
	
	rho0 = 1./(1 + r.^2./8).^2;
	rho0M = diag(rho0);
	g0 = r./(2.*(r.^2./8 + 1));
	g0M = diag(g0);
	zer = zeros(sz);
	
	%M11 = -k^2*ident;
	M12 = -(deriv + rInvM);
	%M13 = -k^2*rho0M;
	M21 = deriv + g0M;
	M22 = ident;
	M23 = rho0M*deriv;
	M31 = ident;
	M32 = zer;
	%M33 = -(deriv2 + rInvM*deriv - k^2*ident);
	
	B = [ident zer zer
		zer zer zer
		zer zer zer];
	
	omg = zeros(numKVals,1);
	eigVects = zeros(3*sz,numKVals);
	
	for i=1:numKVals
		
		k = kVals(i);
		M11 = -k^2*ident;
		M13 = -k^2*rho0M;
		M33 = -(deriv2 + rInvM*deriv - k^2*ident);

		A = [M11 M12 M13
			M21 M22 M23
			M31 M32 M33];

		[V,D] = eig(A,B);
		[D,I] = max(abs(D));
		
		
		eigVects(:,i) = V(:,I(1));
		omg(i) = D(1);
		disp(i)
	
	end
	
	
	
	
	
end



