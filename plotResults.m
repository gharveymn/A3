function plotResults(eigVects,omg,kVals,clrs,bounds,plotVects,figs)
	
	omgR = real(omg);
	omgI = imag(omg);
	
	minr = bounds(1);
	maxr = bounds(2);
	
	[oM,oMInd] = max(omgR);
	
	if(plotVects)
	
		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);		
		
		%Begin plotting
		for i=1:numKVals
			
			if(par.useScatter)
				set(groot,'currentfigure',figs{1})
				scatter3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),5,clrs(i),'.');
				set(groot,'currentfigure',figs{2})
				scatter3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),5,clrs(i),'.');
				set(groot,'currentfigure',figs{3})
				scatter3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),5,clrs(i),'.');
			else
				set(groot,'currentfigure',figs{1})
				plot3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),'Color',clrs(i));
				set(groot,'currentfigure',figs{2})
				plot3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),'Color',clrs(i));
				set(groot,'currentfigure',figs{3})
				plot3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),'Color',clrs(i));
			end
			
		end
		
		%\rho_1
		set(groot,'currentfigure',figs{1})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\rho_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\rho_1)$','interpreter','latex','fontsize',14);
		title('$\rho_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		drawnow;
		
		%s_1
		set(groot,'currentfigure',figs{2})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(s_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(s_1)$','interpreter','latex','fontsize',14);
		title('$s_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		drawnow;

		
		%\phi_1
		set(groot,'currentfigure',figs{3})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\phi_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\phi_1)$','interpreter','latex','fontsize',14);
		title('$\phi_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		%Trouble with boundaries messing up plots, 5std should do the trick
		yl = 5*std(reshape(eigVectsR(2*sz+1:end,:),numKVals*sz,1));
		ylim([-yl,yl]);
		drawnow;
		
		
		%k vs real and imaginary -\omega^2
		set(groot,'currentfigure',figs{4})
		scatter3(kVals,omgR,omgI,[],clrs,'.')
		hold on
		scatter3(kVals(oMInd),oM,omgI(oMInd),[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('$-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;

		%k vs imaginary part of -\omega^2
		set(groot,'currentfigure',figs{5})
		scatter(kVals,omgI,[],'b','.')
		hold on
		scatter(kVals(oMInd),oM,[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
		%k vs real part of -\omega^2
		set(groot,'currentfigure',figs{6})
		scatter(kVals,omgR,[],clrs,'.');
		hold on
		scatter(kVals(oMInd),oM,[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
	else
		
		%Same stuff but without plotting the vectors
				
				%k vs real and imaginary -\omega^2
		set(groot,'currentfigure',figs{4})
		scatter3(kVals,omgR,omgI,[],clrs,'.')
		hold on
		scatter3(kVals(oMInd),oM,omgI(oMInd),[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('$-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;

		%k vs imaginary part of -\omega^2
		set(groot,'currentfigure',figs{5})
		scatter(kVals,omgI,[],'b','.')
		hold on
		scatter(kVals(oMInd),oM,[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
		
		%k vs real part of -\omega^2
		set(groot,'currentfigure',figs{6})
		scatter(kVals,omgR,[],clrs,'.');
		hold on
		scatter(kVals(oMInd),oM,[],clrs(oMInd),'.');
		hold off
		xlabel('$k$','interpreter','latex','fontsize',14);
		ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
		title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
		drawnow;
	
	end
end

