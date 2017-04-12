function P2T3
	
	par = getPar;
	
	figs = par.figs;
	
	if(par.MP)
		mp.Digits(34);
	
		numKVals = 10;
		kMin = mp('0');
		kMax = mp('1');
		kVals = mp(linspace(kMin,kMax,numKVals))';
		%kVals = logis(kVals,kMin,kMax);
		h = mp('0.1');
		%interval = mp('[0.001,0.021]');
		domain = mp('[0.001,15]');
		numEigs = 1;
		[omg,eigVects,r] = getDataMP(numKVals,kVals,h,domain,numEigs);
	else
		numKVals = par.numKVals;
		kMin = par.kMin;
		kMax = par.kMax;
		kVals = linspace(kMin,kMax,numKVals)';
		%kVals = logis(kVals,kMin,kMax);
		%interval = mp('[0.001,0.021]');
		logisParams = [{par.curveParam},par.domain,par.transformeddomain];
		numEigs = 1;
		[omg,eigVects,r] = getData(numKVals,kVals,logisParams,numEigs,par.numPoints);
	end
	
	%disp(r(1:100))
	%disp(omg)
	kPlots = kron(kVals,ones(numEigs,1));
	
	szEV = size(eigVects,2);
	sz = size(r,1);
	
	
	if(par.plotImg)
	
		omgR = real(omg);
		omgI = imag(omg);

		%imMin = min(omgI);
		%imMax = max(omgI);
		
		[oM,oMInd] = max(omgR);
		disp(['Max -omega^2 found at k = ' num2str(kVals(oMInd))])
		disp(['Largest value for -omega^2: ' num2str(oM)])

		normalize = kMax-kMin;
		if(~normalize)
			normalize = 1;
		end

		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);
		
		cRs = abs((kMax - kVals)/normalize);

		set(groot,'currentfigure',figs{1})
		clf
		hold on
		for i=1:szEV
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			
			if(par.useScatter)
				scatter3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),5,clr,'.');
			else
				plot3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),'Color',clr);
			end
			
		end
		%plot(r,eigVects(1:sz,15))
		%hold on
		%plot(r,eigVects(1:sz,5))
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\rho_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\rho_1)$','interpreter','latex','fontsize',14);
		title('$\rho_1$','interpreter','latex','fontsize',14);
		xlim([0,max(real(r))])
		drawnow;

		set(groot,'currentfigure',figs{2})
		clf
		hold on
		for i=1:szEV
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			
			if(par.useScatter)
				scatter3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),5,clr,'.');
			else
				plot3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),'Color',clr);
			end
			
		end
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(s_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(s_1)$','interpreter','latex','fontsize',14);
		title('$s_1$','interpreter','latex','fontsize',14);
		xlim([0,max(real(r))])
		%plot(r,eigVects(sz+1:2*sz,15))
		drawnow;

		set(groot,'currentfigure',figs{3})
		clf
		hold on
		for i=1:szEV
			if(i == oMInd)
				clr = [0,1,0];	
			else
				clr = [1-cRs(i),0,cRs(i)];
			end
			if(par.useScatter)
				scatter3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),5,clr,'.');
			else
				plot3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),'Color',clr);
			end
		end
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\phi_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\phi_1)$','interpreter','latex','fontsize',14);
		title('$\phi_1$','interpreter','latex','fontsize',14);
		xlim([0,max(real(r))])
		%plot(r,eigVects(2*sz+1:end,15))
		drawnow;
		
		set(groot,'currentfigure',figs{4})
		clf
		hold on
		%plot(r(1:size(negomgsq)),negomgsq);
		cR = kron(abs((kMax - kVals)./normalize),ones(numEigs,1));
		clr = [1-cR,zeros(size(cR,1),1),cR];
		clr((oMInd-1)*numEigs + 1:oMInd*numEigs,:) = repmat([0,1,0],numEigs,1);
		scatter3(kPlots,omgR,omgI,[],clr,'.')
		scatter3(kPlots(oMInd),oM,omgI(oMInd),[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('$-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;

		set(groot,'currentfigure',figs{5})
		clf
		hold on
		scatter(kPlots,omgI,[],'b','.')
		scatter(kPlots(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
		set(groot,'currentfigure',figs{6})
		clf
		hold on
		scatter(kPlots,omgR,[],clr,'.');
		scatter(kPlots(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
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

