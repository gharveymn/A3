function [M12,M23,SW,dg1du,d2g1du2,uInvM,ident,rho0M]=BuildMatrices(r,u,sz,rho0,g0,logisParams,dtransform)
	
	ident = speye(sz);
	diagn = (1:sz)';
	uInvM = sparse(diagn,diagn,1./u);
	
	i = [(2:sz-1)';(2:sz-1)'];
	j = [(3:sz)';(1:sz-2)'];
	v = [ones(sz-2,1);-ones(sz-2,1)];
	
	%Derivatives:
	shftforward2 = circshift(r,-2);
	preDif = shftforward2(1:end-2)-r(1:end-2);
	difr = [preDif;preDif];
	
	ddr = sparse(i,j,v./difr,sz,sz);
	ddr(1,1) = -1/(r(2)-r(1));
	ddr(1,2) = -ddr(1,1);
	ddr(sz,sz) = 1/(r(end)-r(end-1));
	ddr(sz,sz-1) = -ddr(sz,sz);
	
	%Asymmetric second derivative operator
	%	MatRow <- [2/(h_1*(h_1+h_2)) , -2*f(x)/(h_1*h_2) , 2/(h_2*(h_1+h_2))
	
	i = [(2:sz)';(1:sz)';(1:sz-1)'];
	j = [(1:sz-1)';(1:sz)';(2:sz)'];
	
	shftback = circshift(r,1);
	shftforward = circshift(r,-1);
	h1vec = r(2:end-1)-shftback(2:end-1);		%for easier creation of diags
	h2vec = shftforward(2:end-1)-r(2:end-1);
	h1veclast = r(end)-shftback(end);
	h2vecfirst = shftforward(1)-r(1);
	
	diagleft = sparse((1:sz-2)',1,2./(h1vec.*(h1vec+h2vec)),sz-1,1);
	diagleft(end) = 1./(h1veclast.^2);
	
	diagmiddle = sparse((2:sz-1)',1,-2./(h1vec.*h2vec),sz,1);
	diagmiddle(1) = -2./(h2vec(1).^2);
	diagmiddle(end) = -2./(h1vec(end).^2);
	
	diagright = sparse((2:sz-1)',1,2./(h2vec.*(h1vec+h2vec)),sz-1,1);
	diagright(1) = 1./(h2vecfirst.^2);
	
	difr = [diagleft;diagmiddle;diagright];
	
	d2fdr2 = sparse(i,j,difr,sz,sz);
	
	
	drdu = sparse(diagn,diagn,...
				dtransform(1,u,logisParams{:})...
				);
				
	d2rdu2 = sparse(diagn,diagn,...
				dtransform(2,u,logisParams{:})...
				);
	
	
	%plot(r,diag(drdu))
	%plot(r,diag(d2rdu2))
	
	dfdu = drdu*ddr;
	
	rho0M = sparse(diagn,diagn,rho0);
	g0M = sparse(diagn,diagn,g0);
	zer = sparse(sz,sz);
	
	dg1du = dfdu;
	%dg1du(1,:) = 0;
	dg1du(end,:) = 0;
	%dg1du(:,1) = 0;
	%dg1du(:,end) = 0;
	
	dg1dr = ddr;
	%dg1dr(1,:) = 0;
	%dg1dr(end,:) = 0;
	%dg1dr(:,1) = 0;
	%dg1dr(:,end) = 0;
	
	d2g1dr2 = d2fdr2;
	%d2g1dr2(1,:) = 0;
	%d2g1dr2(end,:) = 0;
	%d2g1dr2(:,1) = 0;
	%d2g1dr2(:,end) = 0;
	
	d2g1du2 = (drdu^2)*d2g1dr2 + d2rdu2*dg1dr;
	%d2g1du2(1,:) = 0;
	%d2g1du2(end,:) = 0;
	%d2g1du2(:,1) = 0;
	%d2g1du2(:,end) = 0;
	
	%M11 = -k^2*ident;
	M12 = -(dfdu + uInvM);
	%M13 = -k^2*rho0M;
	%M21 = deriv + g0M;
	%M22 = ident;
	M23 = dg1du*rho0M;	%I know this isn't right, but this is the only way I can get the code to behave.
	%M31 = ident;
	%M32 = zer;
	%M33 = -(deriv2 + rInvM*deriv - k^2*ident);
	
	%Southwest corner of the matrix, [M21 M22; M31 M32].
	SW = [dfdu+g0M, ident
		ident,zer];
	
end