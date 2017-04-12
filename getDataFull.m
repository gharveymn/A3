function [omg,eigVects,u,conds]=getDataFull(numKVals,kVals,domain,transformeddomain,numEigs,numPoints)
		
	%Build Matrices
	r = linspace(domain{:},numPoints)';
	%bilogis is centered at midpoints of domain and range, this finds half of domain and range
	logitParams = {10,domain{:},transformeddomain{:}};
	u = halfbilogit(r,logitParams{:});
	sz = size(u,1);
	
	rho0 = 1./(1 + r.^2./8).^2;
	g0 = r./(2.*(r.^2./8 + 1));
	
	[M12,M23,SW,dfdr,d2fdr2,dg1dr,d2g1dr2,rInvM,uInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0,logitParams);
	B = diag([ones(1,sz),ones(1,2*sz)]);
	
	omg = zeros(numEigs*numKVals,1);
	eigVects = zeros(3*sz,numEigs*numKVals);
	
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
		
		[V,D] = eig(A,B,'vector');
		
		DF = D(~imag(D) & isfinite(D));
		VF = V(:,~imag(D) & isfinite(D));
		[minDF,minDFI] = min(DF);
		omg(i) = minDF;
		eigVects(:,i) = VF(:,minDFI);
		disp(i)
	
	end
	
	
	
	
end

function [M12,M23,SW,dfdr,d2fdr2,dg1dr,d2g1dr2,rInvM,uInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0,logitParams)
	
	ident = eye(sz);
	rInvM = diag(1./r);
	uInvM = diag(1./u);
	
	%Derivatives:
	shft = circshift(u,-2);
	preDif = shft(1:end-2)-u(1:end-2);
	
	ddu = diag(preDif,1) - diag(preDif,-1);
	ddu(1,1) = -1/(u(2)-u(1));
	ddu(1,2) = -ddu(1,1);
	ddu(sz,sz) = 1/(u(end)-u(end-1));
	ddu(sz,sz-1) = -ddu(sz,sz);
	
	d2du2 = ddu^2;
	
	dudr = diag(dhalfbilogit(r,1,logitParams{:}));
	d2udr2 = diag(dhalfbilogit(r,2,logitParams{:}));
	%d2dudr(1,:) = zeros(1,sz);
	%d2dudr(end,:) = d2dudr(1,:);
	
	dfdr = dudr*ddu;
	d2fdr2 = (dudr.^2)*d2du2 + d2udr2*ddu;
	%d2fdr2 = d2udr2*ddu;
	
	rho0M = diag(rho0);
	%g0 = 4.*r./(8 + r.^2);
	g0M = diag(g0);
	zer = zeros(sz);
	
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
	M12 = -(dfdr + rInvM);
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




