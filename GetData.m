function [omg,eigVects,r]=GetData(numKVals,kVals,logisParams,numPoints,transform,dtransform,plotPack)
	
	
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
	
	%plot(r,u)
	
	rho0 = 1./(1 + u.^2./8).^2;
	g0 = u./(2.*(u.^2./8 + 1));
	
	[M12,M23,SW,dg1du,d2g1du2,uInvM,ident,rho0M]=BuildMatrices(r,u,sz,rho0,g0,logisParams,dtransform);
	B = sparse(1:sz,1:sz,ones(1,sz),3*sz,3*sz);
	
	omg = zeros(numKVals,1);
	
	if(plotVects)
		eigVects = zeros(3*sz,numKVals);
	else
		eigVects = [];
	end
	
	v0 = [rho0;r*.03;-log(rho0)];
	v0(sz) = 0; %Apply rho1 boundary cond
	v0(sz+1) = 0; %Apply sr boundary cond
	
	%figure
	%hold on
	for i=1:numKVals
		
		k = kVals(i);
		M11 = -k^2*ident;
		M13 = -k^2*rho0M;
		M33 = -d2g1du2 - uInvM*dg1du + k^2*ident;
		%M33 = -d2g1dr2 + k^2*ident;

		A = [M11 M12 M13
			SW,[M23;M33]];
		A(2*sz+1,:) = 0;
		A(2*sz+1,2*sz+1) = 1;
		A(2*sz+1,2*sz+2) = -1;
		
		opts = struct('isreal',1,'p',50,'maxit',200,'disp',0,'v0',v0);
		
		try
			[V,D] = eigs(A,B,3*sz,'sm',opts);
			
			D = diag(D);
			DF = D(~imag(D) & isfinite(D));
			[maxDF,maxDFI] = max(DF);
			omg(i) = maxDF;
		
			if(livePlot)
				if(plotVects)
					VF = V(:,~imag(D) & isfinite(D));
					v = VF(:,maxDFI);
					eigVects(:,i) = v;
					PlotResults(r,v,omg(1:i),kVals(1:i),plotPack{1:end-2},plotPack{end-1}(1:i,:))
				else
					PlotResults(r,[],omg(1:i),kVals(1:i),plotPack{1:end-2},plotPack{end-1}(1:i,:))
				end

			end
			
			disp([num2str(i) ': k = ' num2str(k)]);
		
		catch
			%warning('There was an error with eigs')
			%disp('caught')
		end
	
	end
	
end




