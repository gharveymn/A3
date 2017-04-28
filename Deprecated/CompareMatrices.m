function [GDA,GDB] = CompareMatrices
	
	par = Parameters;
	logisParams = [{par.curveParam},par.domain,par.transformeddomain];
	numPoints = 10;
	
	k = 0.1;
	
	r = linspace(log(logisParams{2}),log(logisParams{3}),numPoints)';
	
	u = par.T(r,logisParams{:});
	sz = numel(u);
	
	%plot(r,u)
	
	rho0 = 1./(1 + u.^2./8).^2;
	g0 = u./(2.*(u.^2./8 + 1));
	
	[M12,M23,SW,dg1du,d2g1du2,uInvM,ident,rho0M]=BuildMatrices(r,u,sz,rho0,g0,logisParams,par.dT);
	GDB = sparse(1:sz,1:sz,ones(1,sz),3*sz,3*sz);
	
	M11 = -k^2*ident;
	M13 = -k^2*rho0M;
	M33 = -d2g1du2 - uInvM*dg1du + k^2*ident;
	%M33 = -d2g1dr2 + k^2*ident;

	GDA = [M11 M12 M13
		SW,[M23;M33]];
end

