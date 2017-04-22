function PlotResults(r,eigVects,omg,kVals,minr,maxr,plotVects,useScatter,figs,clrs)
		
	omgR = real(omg);
	omgI = imag(omg);
	
	[oM,oMInd] = max(omgR);
	
	oldColor = clrs(oMInd,:);
	clrs(oMInd,:) = [0,1,0];
	
	sz = numel(r);
	
	if(plotVects)
	
		eigVectsR = real(eigVects);
		eigVectsI = imag(eigVects);
		
		
		%Quickfix, messy...	
		%It runs slow if I don't do this
		plotOne = (size(eigVects,2) == 1);
		
		if(plotOne)
			%Access the current color index if we only plot one vector
			to = 1;
		else
			to = numel(kVals);
		end
			
		%Begin plotting
		for i=1:to
			
			if(plotOne)
				clrs(oMInd,:) = oldColor;
				clr = clrs(numel(kVals),:);
				
				%reset for -\omega^2 plotting
				clrs(oMInd,:) = [0,1,0];
			else
				clr = clrs(i,:);
			end
			
			if(useScatter)
				set(groot,'currentfigure',figs{1})
				scatter3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),5,clr,'.');
				set(groot,'currentfigure',figs{2})
				scatter3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),5,clr,'.');
				set(groot,'currentfigure',figs{3})
				scatter3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),5,clr,'.');
			else
				set(groot,'currentfigure',figs{1})
				plot3(r,eigVectsR(1:sz,i),eigVectsI(1:sz,i),'Color',clr);
				set(groot,'currentfigure',figs{2})
				plot3(r,eigVectsR(sz+1:2*sz,i),eigVectsI(sz+1:2*sz,i),'Color',clr);
				set(groot,'currentfigure',figs{3})
				plot3(r,eigVectsR(2*sz+1:end,i),eigVectsI(2*sz+1:end,i),'Color',clr);
			end
			
		end
		
		az = 0;
		el = 90;
		
		%\rho_1
		set(groot,'currentfigure',figs{1})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\rho_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\rho_1)$','interpreter','latex','fontsize',14);
		title('$\rho_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		view(az, el);
		drawnow;
		
		%s_1
		set(groot,'currentfigure',figs{2})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(s_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(s_1)$','interpreter','latex','fontsize',14);
		title('$s_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		view(az, el);
		drawnow;

		
		%\phi_1
		set(groot,'currentfigure',figs{3})
		xlabel('$r$','interpreter','latex','fontsize',14);
		ylabel('Re$(\phi_1)$','interpreter','latex','fontsize',14);
		zlabel('Im$(\phi_1)$','interpreter','latex','fontsize',14);
		title('$\phi_1$','interpreter','latex','fontsize',14);
		xlim([minr,maxr])
		%Trouble with boundaries messing up plots, 5std should do the trick
		%yl = 5*std(reshape(eigVectsR(2*sz+1:end,1:),numel(kVals)*sz,1));
		%ylim([-yl,yl]);
		view(az, el);
		drawnow;
	end
	
	%Plot the -\omega^2 values wrt k
	
	%k vs real and imaginary -\omega^2
	set(groot,'currentfigure',figs{4})
	clf
	scatter3(kVals,omgR,omgI,[],clrs,'.')
	hold on
	scatter3(kVals(oMInd),oM,omgI(oMInd),[],clrs(oMInd,:),'.');
	hold off
	xlabel('$k$','interpreter','latex','fontsize',14);
	ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
	zlabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
	title('$-\omega^2$','interpreter','latex','fontsize',14);
	drawnow;

	%k vs imaginary part of -\omega^2
	set(groot,'currentfigure',figs{5})
	clf
	scatter(kVals,omgI,[],'b','.')
	hold on
	scatter(kVals(oMInd),omgI(oMInd),[],clrs(oMInd,:),'.');
	hold off
	xlabel('$k$','interpreter','latex','fontsize',14);
	ylabel('Im$(-\omega^2)$','interpreter','latex','fontsize',14);
	title('Imaginary part of $-\omega^2$','interpreter','latex','fontsize',14);
	drawnow;

	%k vs real part of -\omega^2
	set(groot,'currentfigure',figs{6})
	clf
	scatter(kVals,omgR,[],clrs,'.');
	hold on
	scatter(kVals(oMInd),oM,[],clrs(oMInd,:),'.');
	hold off
	xlabel('$k$','interpreter','latex','fontsize',14);
	ylabel('Re$(-\omega^2)$','interpreter','latex','fontsize',14);
	title('Real part of $-\omega^2$','interpreter','latex','fontsize',14);
	drawnow;
	
end

