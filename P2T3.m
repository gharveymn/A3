function figs=P2T3
	
	%Get parameters
	par = getPar;
	
	numKVals = par.numKVals;
	kMin = par.kMin;
	kMax = par.kMax;
	kVals = linspace(kMin,kMax,numKVals)';
	logisParams = [{par.curveParam},par.domain,par.transformeddomain];
	
	if(par.plotVects)
		figs = {figure(1),figure(2),figure(3),figure(4),figure(5),figure(6)};
	else
		figs = {figure(1),figure(2),figure(3)};
	end
	drawnow
	
	%endblock
	
	%*Data collection*
	[omg,eigVects,r] = getData(numKVals,kVals,logisParams,par.numPoints,par.plotVects);
	
	
	szEV = size(eigVects,2);
	sz = size(r,1);
	
	omgR = real(omg);
	omgI = imag(omg);

	[oM,oMInd] = max(omgR);
	disp(['Max -omega^2 found at k = ' num2str(kVals(oMInd))])
	disp(['Largest value for -omega^2: ' num2str(oM)])


	%Create color map
	normalize = kMax-kMin;
	if(~normalize)
		normalize = 1;
	end
	cRs = abs((kMax - kVals)./normalize);
	
	
	if(par.plotVects)
	
		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);		
		
		%Begin plotting
		
		%\rho_1
		set(groot,'currentfigure',figs{1})
		clf
		hold on
		for i=1:szEV
			
			%If this is the max -\omega^2 value, plot green, otherwise continue with R-B gradient
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
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\rho_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\rho_1)$','interpreter','latex','fontsize',14);
		title('$\rho_1$','interpreter','latex','fontsize',14);
		xlim([0,max(real(r))])
		drawnow;

		
		%s_1
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
		drawnow;

		
		%\phi_1
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
		
		%Trouble with boundaries messing up plots, 5std should do the trick
		yl = 5*std(reshape(eigVectsR(2*sz+1:end,:),szEV*sz,1));
		ylim([-yl,yl]);
		drawnow;
		
		
		%k vs real and imaginary -\omega^2
		set(groot,'currentfigure',figs{4})
		clf
		hold on
		
		%Color mapping
		clr = [1-cRs,zeros(size(cRs,1),1),cRs];
		scatter3(kVals,omgR,omgI,[],clr,'.')
		scatter3(kVals(oMInd),oM,omgI(oMInd),[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('$-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;

		%k vs imaginary part of -\omega^2
		set(groot,'currentfigure',figs{5})
		clf
		hold on
		scatter(kVals,omgI,[],'b','.')
		scatter(kVals(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
		%k vs real part of -\omega^2
		set(groot,'currentfigure',figs{6})
		clf
		hold on
		scatter(kVals,omgR,[],clr,'.');
		scatter(kVals(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
	else
		
		%Same stuff but without plotting the vectors
				
		%k vs real and imaginary -\omega^2
		set(groot,'currentfigure',figs{1})
		clf
		hold on
		
		%Color mapping
		clr = [1-cRs,zeros(size(cRs,1),1),cRs];
		scatter3(kVals,omgR,omgI,[],clr,'.')
		scatter3(kVals(oMInd),oM,omgI(oMInd),[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('$-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;

		%k vs imaginary part of -\omega^2
		set(groot,'currentfigure',figs{2})
		clf
		hold on
		scatter(kVals,omgI,[],'b','.')
		scatter(kVals(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
		%k vs real part of -\omega^2
		set(groot,'currentfigure',figs{3})
		clf
		hold on
		scatter(kVals,omgR,[],clr,'.');
		scatter(kVals(oMInd),oM,[],[0,1,0],'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
	
	end
	
	
end

