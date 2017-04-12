function [omg,eigVects,u,conds]=getData(numKVals,kVals,h,interval,numEigs,numPoints)
		
	%Build Matrices
	r = linspace(interval(1),interval(2),numPoints)';
	u = log(r);
	plot(r,u)
	drawnow
	%r
	%r = logis(r,interval(1),interval(2));
	sz = size(u,1);
	
	rho0 = 1./(1 + r.^2./8).^2;
	g0 = r./(2.*(r.^2./8 + 1));
	
	[M12,M23,SW,dfdr,d2fdr2,dg1dr,d2g1dr2,uInvM,rInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0);
	B = sparse(1:sz,1:sz,ones(1,sz),3*sz,3*sz);
	
	omg = zeros(numEigs*numKVals,1);
	eigVects = zeros(3*sz,numEigs*numKVals);
	
	v0 = [rho0;u*.03;-log(rho0)];
	v0(sz) = 0; %Apply rho1 boundary cond
	v0(sz+1) = 0; %Apply sr boundary cond
	
	conds = zeros(size(kVals,1),1);
	%figure
	%hold on
	for i=1:numKVals
		
		k = kVals(i);
		M11 = -k^2*ident;
		M13 = -k^2*rho0M;
		M33 = -(d2g1dr2 + rInvM*dg1dr - k^2*ident);

		A = [M11 M12 M13
			SW,[M23;M33]];
		
		
		opts = struct('isreal',1,'v0',v0);
		
		[V,D] = eigs(A,B,1,'sm',opts);
		D = diag(D);
		omg(i) = D(1);
		eigVects(:,i) = V(:,1);
		disp(i)
	
	end
	
	
	
	
end

function [M12,M23,SW,dfdr,d2fdr2,dg1dr,d2g1dr2,uInvM,rInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0)
	
	ident = speye(sz);
	diagn = (1:sz)';
	uInvM = sparse(diagn,diagn,1./u);
	rInvM = sparse(diagn,diagn,1./r);
	
	i = [(2:sz-1)';(2:sz-1)'];
	j = [(3:sz)';(1:sz-2)'];
	v = [ones(sz-2,1);-ones(sz-2,1)];
	
	
	%Derivatives:
	shft = circshift(u,-2);
	preDif = shft(1:end-2)-u(1:end-2);
	difu = [preDif;preDif];
	
	ddu = sparse(i,j,v./difu,sz,sz);
	ddu(1,1) = -1/(u(2)-u(1));
	ddu(1,2) = -ddu(1,1);
	ddu(sz,sz) = 1/(u(end)-u(end-1));
	ddu(sz,sz-1) = -ddu(sz,sz);
	
	d2du2 = ddu^2;
	
	dudr = uInvM;
	d2udr2 = -uInvM.^2;
	%d2dudr(1,:) = zeros(1,sz);
	%d2dudr(end,:) = d2dudr(1,:);
	
	dfdr = dudr*ddu;
	d2fdr2 = (dudr.^2)*d2du2 + d2udr2*ddu;
	%d2fdr2 = d2udr2*ddu;
	
	rho0M = sparse(diagn,diagn,rho0);
	%g0 = 4.*r./(8 + r.^2);
	g0M = sparse(diagn,diagn,g0);
	zer = sparse(sz,sz);
	
	dg1dr = dfdr;
	dg1dr(1,:) = zeros(1,sz);
	dg1dr(end,:) = dg1dr(1,:);
	
	dg1du = ddu;
	dg1du(1,:) = zeros(1,sz);
	dg1du(end,:) = dg1du(1,:);
	
	d2g1du2 = dg1du^2;
	
	d2g1dr2 = (dudr.^2)*d2g1du2 + d2udr2*dg1du;
	%d2g1dr2 = d2udr2*dg1du;
	
	%M11 = -k^2*ident;
	M12 = -(dfdr + uInvM);
	%M13 = -k^2*rho0M;
	%M21 = deriv + g0M;
	%M22 = ident;
	M23 = rho0M*dg1dr;
	%M31 = ident;
	%M32 = zer;
	%M33 = -(deriv2 + rInvM*deriv - k^2*ident);
	
	SW = [dfdr+g0M, ident
		ident,zer];
	
end




