function [omg,eigVects,r,conds]=getData(numKVals,kVals,h,interval,numEigs)
		
	%Build Matrices
	r = (interval(1):h:interval(2))';
	%r
	%r = logis(r,interval(1),interval(2));
	sz = size(r,1);
	
	rho0 = 1./(1 + r.^2./8).^2;
	g0 = r./(2.*(r.^2./8 + 1));
	
	[M12,M23,SW,deriv,deriv2,g1Deriv,rInvM,ident,rho0M]=buildMatrices(r,sz,rho0,g0);
	B = sparse(1:sz,1:sz,ones(1,sz),3*sz,3*sz);
	
	omg = zeros(numEigs*numKVals,1);
	eigVects = zeros(3*sz,numEigs*numKVals);
	
	v0 = [rho0;r*.03;-log(rho0)];
	v0(sz) = 0; %Apply rho1 boundary cond
	v0(sz+1) = 0; %Apply sr boundary cond
	
	conds = zeros(size(kVals,1),1);
	%figure
	%hold on
	for i=1:numKVals
		
		k = kVals(i);
		M11 = -k^2*ident;
		M13 = -k^2*rho0M;
		M33 = -(g1Deriv^2 + rInvM*g1Deriv - k^2*ident);

		A = [M11 M12 M13
			SW,[M23;M33]];
		
		opts = struct('isreal',1,'v0',v0);
		
		
		[V,D] = eigs(A,B,5,'sm',opts);
		D = diag(D);
		omg(i) = D(5);
		eigVects(:,i) = V(:,5);
		%DF = D(isfinite(D));
		%VF = V(:,isfinite(D));
		%[maxDF,maxDFI] = max(real(DF));
		%omg(i) = DF(maxDFI);
		%eigVects(:,i) = VF(:,maxDFI);
		
		
		%{
		for j = 1:20
			[V,D] = eigs(A,B,j,'sm',opts);
			disp(['(' num2str(i) ',' num2str(j) ')'])
			disp(D(end))
			scatter(j,real(D(end)),10,[1-1/i,0,1/i],'.')
			drawnow
		end
		%}
		%{
		[V,D] = eigs(A,B,numEigs,'sm',opts);
		j = 2;
		
		while(imag(D(end)) ~= 0)
			[V,D] = eigs(A,B,j,'sm',opts);
			disp(['(' num2str(i) ',' num2str(j) ')'])
			disp(D(end))
			j = j+1;
		end
		%}
		%omg(i) = D(end);
		%eigVects(:,i) = V(:,end);
		disp(i)
	
	end
	
	
	
end

function [M12,M23,SW,deriv,deriv2,g1Deriv,rInvM,ident,rho0M]=buildMatrices(r,sz,rho0,g0)
	
	%TODO: Fix derivative matrices 
	shft = circshift(r,-2);
	preDif = shft(1:end-2)-r(1:end-2);
	difR = [preDif;preDif];
	
	rInv = 1./r;
	diagn = (1:sz)';
	rInvM = sparse(diagn,diagn,rInv);
	%rInvM = sparse(diag(rInv));
	ident = speye(sz);
	
	i = [(2:sz-1)';(2:sz-1)'];
	j = [(3:sz)';(1:sz-2)'];
	v = [ones(sz-2,1);-ones(sz-2,1)];
	
	deriv = sparse(i,j,v./difR,sz,sz);
	%deriv = sparse(toeplitz([0;-1;zeros(sz-2,1)],[0,1,zeros(1,sz-2)]));
	%deriv = deriv./(2*h);
	deriv(1,1) = -1/(r(2)-r(1));
	deriv(1,2) = -deriv(1,1);
	deriv(sz,sz) = 1/(r(end)-r(end-1));
	deriv(sz,sz-1) = -deriv(sz,sz);
	
	deriv2 = deriv^2;
	deriv2(1,:) = zeros(1,sz);
	deriv2(end,:) = deriv2(1,:);
	
	rho0M = sparse(diagn,diagn,rho0);
	%g0 = 4.*r./(8 + r.^2);
	g0M = sparse(diagn,diagn,g0);
	zer = sparse(sz,sz);
	
	g1Deriv = deriv;
	g1Deriv(1,:) = zeros(1,sz);
	g1Deriv(end,:) = g1Deriv(1,:);
	
	%M11 = -k^2*ident;
	M12 = -(deriv + rInvM);
	%M13 = -k^2*rho0M;
	%M21 = deriv + g0M;
	%M22 = ident;
	M23 = rho0M*g1Deriv;
	%M31 = ident;
	%M32 = zer;
	%M33 = -(deriv2 + rInvM*deriv - k^2*ident);
	
	SW = [deriv+g0M, ident
		ident,zer];
	
end




