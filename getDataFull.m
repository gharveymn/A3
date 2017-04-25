function [omg,eigVects,r]=GetDataFull(numKVals,kVals,logisParams,numPoints,transform,dtransform,plotPack)
		
	
	plotVects = plotPack{3};
	livePlot = plotPack{7};
	
	r = linspace(logisParams{2},logisParams{3},numPoints)';
	
	
	badzero = find(~r);
	if(~isempty(badzero))
		r(badzero) = r(badzero+1)/2;
		if(max(r) == -min(r))
			disp('Warning: Rectified zero in domain; domain is not symmetric.')
		else
			disp('Warning: Rectified zero in domain')
		end
	end
	
	u = transform(r,logisParams{:});
	sz = numel(u);
	
	rho0 = 1./(1 + r.^2./8).^2;
	g0 = r./(2.*(r.^2./8 + 1));
	
	[M12,M23,SW,dg1dr,d2g1dr2,rInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0,logisParams,dtransform);
	B = diag([ones(1,sz),ones(1,2*sz)]);
	
	omg = zeros(numKVals,1);
	if(plotVects)
		eigVects = zeros(3*sz,numKVals);
	else
		eigVects = [];
	end
	
	%figure
	%hold on
	for i=1:numKVals
		
		k = kVals(i);
		M11 = -k^2*ident;
		M13 = -k^2*rho0M;
		M33 = -(d2g1dr2 + rInvM*dg1dr - k^2*ident);

		A = [M11 M12 M13
			SW,[M23;M33]];
		
		
		try
			[V,D] = eig(A,B,'vector');
		
			DF = D(~imag(D) & isfinite(D));
			[minDF,minDFI] = min(DF);
			omg(i) = minDF;
		
			if(livePlot)
				if(plotVects)
					VF = V(:,~imag(D) & isfinite(D));
					v = VF(:,minDFI);
					eigVects(:,i) = v;
					PlotResults(r,v,omg(1:i),kVals(1:i),plotPack{1:end-2},plotPack{end-1}(1:i,:))
				else
					PlotResults(r,[],omg(1:i),kVals(1:i),plotPack{1:end-2},plotPack{end-1}(1:i,:))
				end
			end
			
			disp([num2str(i) ': k = ' num2str(k)]);
		
		catch
			warning('There was an error with eig')
		end
	
	end
	
	
	
	
end

function [M12,M23,SW,dg1dr,d2g1dr2,rInvM,ident,rho0M]=buildMatrices(r,u,sz,rho0,g0,logisParams,dtransform)
	
	ident = eye(sz);
	rInvM = diag(1./r);
	
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
	
	ddu = full(ddu);
	
	d2du2 = ddu^2;
	
	dudr = diag(dtransform(1,r,logisParams{:}));
	d2udr2 = diag(dtransform(2,r,logisParams{:}));
	%d2dudr(1,:) = zeros(1,sz);
	%d2dudr(end,:) = d2dudr(1,:);
	
	dfdr = dudr*ddu;
	d2fdr2 = (dudr^2)*d2du2 + d2udr2*ddu;
	%d2fdr2 = d2udr2*ddu;
	
	rho0M = diag(rho0);
	%g0 = 4.*r./(8 + r.^2);
	g0M = diag(g0);
	zer = zeros(sz);
	
	dg1du = ddu;
	%dg1du(1,:) = 0;
	%dg1du(end,:) = 0;
	%dg1du(:,1) = 0;
	%dg1du(:,end) = 0;
	
	dg1dr = dudr*dg1du;
	dg1dr(1,:) = 0;
	dg1dr(end,:) = 0;
	dg1dr(:,1) = 0;
	dg1dr(:,end) = 0;
	
	d2g1du2 = dg1du^2;
	%d2g1du2(1,:) = 0;
	%d2g1du2(end,:) = 0;
	%d2g1du2(:,1) = 0;
	%d2g1du2(:,end) = 0;
	
	d2g1dr2 = (dudr^2)*d2g1du2 + d2udr2*dg1du;
	d2g1dr2(1,:) = 0;
	d2g1dr2(end,:) = 0;
	%d2g1dr2(:,1) = 0;
	%d2g1dr2(:,end) = 0;
	
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




