function P2T3
	
	par = getPar;
	
	numKVals = 1000;
	kMin = .15;%0.06
	kMax = .55;
	kVals = linspace(kMin,kMax,numKVals)';
	%kVals = logis(kVals,kMin,kMax);
	h = 0.001;
	interval = [0.001,20];
	numEigs = 1;
	
	[omg,eigVects,r] = getData(numKVals,kVals,h,interval,numEigs);
	
	kPlots = kron(kVals,ones(numEigs,1));
	
	szEV = size(eigVects,2);
	sz = size(r,1);
	
	if(par.plotImg)
	
		omgR = real(omg);
		omgI = imag(omg);

		imMin = min(omgI);
		imMax = max(omgI);

		normalize = kMax;
		if(~normalize)
			normalize = 1;
		end

		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);
		

		figure(1)
		clf
		hold on
		for i=1:szEV
			%cR = abs(imMin - omgI(i))/normalize;
			cR = abs((kMax - kVals(i))/normalize);
			plot3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),'Color',[1-cR,0,cR]);
		end
		%plot(r,eigVects(1:sz,15))
		%hold on
		%plot(r,eigVects(1:sz,5))
		xlabel('r');
		ylabel('re(p)');
		zlabel('im(p)');
		drawnow;

		figure(2)
		clf
		hold on
		for i=1:szEV
			%cR = abs(imMin - omgI(i))/normalize;
			cR = abs((kMax - kVals(i))/normalize);
			plot3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),'Color',[1-cR,0,cR]);
		end
		xlabel('r');
		ylabel('re(s)');
		zlabel('im(s)');
		%plot(r,eigVects(sz+1:2*sz,15))
		drawnow;

		figure(3)
		clf
		hold on
		for i=1:szEV
			%cR = abs(imMin - omgI(i))/normalize;
			cR = abs((kMax - kVals(i))/normalize);
			plot3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),'Color',[1-cR,zeros(size(cR,1),1),cR]);
		end
		xlabel('r');
		ylabel('re(phi)');
		zlabel('im(phi)');
		%plot(r,eigVects(2*sz+1:end,15))
		drawnow;
		
		figure(4)
		clf
		%plot(r(1:size(negomgsq)),negomgsq);
		cR = kron(abs((kMax - kVals)./normalize),ones(numEigs,1));
		scatter3(kPlots,omgR,omgI,[],[1-cR,zeros(size(cR,1),1),cR],'.')
		%plot(kPlots,abs(omg),'.','Color','r')
		xlabel('k');
		ylabel('re(o)');
		zlabel('im(o)');
		drawnow;
		
		figure(5)
		clf
		scatter(kPlots,omgR,[],[1-cR,zeros(size(cR,1),1),cR],'.');
		xlabel('k');
		makeylabel('-o^2');
		drawnow;

		figure(6)
		clf
		plot(kPlots,omgI,'.','Color','b')
		xlabel('k');
		makeylabel('im(o)');
		drawnow;
	
	else

		figure(1)
		clf
		for i=1:szEV
			plot(r,eigVects(1:sz,i),'.','Color','r');
		end
		xlabel('r');
		ylabel('rho1');
		drawnow;

		figure(2)
		clf
		for i=1:szEV
			plot(r,eigVects(sz+1:2*sz,i),'.','Color','g');
		end
		xlabel('r');
		xlabel('s1');
		drawnow;

		figure(3)
		clf
		for i=1:szEV
			plot(r,eigVects(2*sz+1:end,i),'.','Color','b');
		end
		xlabel('r');
		ylabel('phi1');
		drawnow;
		
		figure(4)
		clf
		for i=1:szEV
			plot(kPlots,omg,'.','Color','k');
		end
		xlabel('k');
		ylabel('-o^2');
		drawnow;
	
	end
	
	
end

