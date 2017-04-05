function P2T3
	
	par = getPar;
	
	if(par.MP)
		mp.Digits(34);
	
		numKVals = 10;
		kMin = mp('0');
		kMax = mp('4');
		kVals = mp(linspace(kMin,kMax,numKVals))';
		%kVals = logis(kVals,kMin,kMax);
		h = mp('0.1');
		%interval = mp('[0.001,0.021]');
		interval = mp('[0.001,15]');
		numEigs = 1;
		[omg,eigVects,r] = getDataMP(numKVals,kVals,h,interval,numEigs);
	else
		numKVals = 100;
		kMin = 0;
		kMax = 4;
		kVals = linspace(kMin,kMax,numKVals)';
		%kVals = logis(kVals,kMin,kMax);
		h = 0.01;
		%interval = mp('[0.001,0.021]');
		interval = [0.001,15];
		numEigs = 1;
		[omg,eigVects,r] = getData(numKVals,kVals,h,interval,numEigs);
	end
	
	disp(omg)
	kPlots = kron(kVals,ones(numEigs,1));
	
	szEV = size(eigVects,2);
	sz = size(r,1);
	
	if(par.plotImg)
	
		omgR = real(omg);
		omgI = imag(omg);

		imMin = min(omgI);
		imMax = max(omgI);
		
		[oM,oMInd] = max(omgR);

		normalize = kMax;
		if(~normalize)
			normalize = 1;
		end

		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);
		
		cRs = abs((kMax - kVals)/normalize);

		figure(1)
		clf
		hold on
		for i=1:szEV
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			scatter3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),5,clr,'.');
			%plot3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),'Color',clr);
			
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
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			scatter3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),5,clr,'.');
			%plot3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),'Color',clr);
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
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			scatter3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),5,clr,'.');
			%plot3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),'Color',clr);
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
		clr = [1-cR,zeros(size(cR,1),1),cR];
		clr((oMInd-1)*numEigs + 1:oMInd*numEigs,:) = repmat([0,1,0],numEigs,1);
		scatter3(kPlots,omgR,omgI,[],clr,'.')
		%plot(kPlots,abs(omg),'.','Color','r')
		xlabel('k');
		ylabel('re(o)');
		zlabel('im(o)');
		drawnow;
		
		figure(5)
		clf
		scatter(kPlots,omgR,[],clr,'.');
		xlabel('k');
		makeylabel('-o^2');
		drawnow;

		figure(6)
		clf
		plot(kPlots,omgI,'.','Color','b')
		xlabel('k');
		makeylabel('im(o)');
		drawnow;
		
		%{
		figure(6)
		clf
		plot(kPlots,conds,'.','Color','k')
		xlabel('k');
		makeylabel('condition number');
		drawnow;
		
		disp(conds)
		%}
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

