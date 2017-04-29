function figs=RunMeSecond
	
	%Get parameters
	
	par = Parameters;
	
	numKVals = par.numKVals;
	kMin = par.kMin;
	kMax = par.kMax;
	kVals = linspace(kMin,kMax,numKVals)';
	logisParams = [{par.curveParam},par.domain,par.transformeddomain];
	numPoints = par.numPoints;
	
	if(par.plotVects)
		
		%Just so that f1,f2,f3 are on top when initializing
		f4 = figure(4);
		f5 = figure(5);
		f6 = figure(6);
		
		f1 = figure(1);
		clf
		hold on
		f2 = figure(2);
		clf
		hold on
		f3 = figure(3);
		clf
		hold on
		
		figs = {f1,f2,f3,f4,f5,f6};
		
	else
		figs = {figure(4),figure(5),figure(6)};
	end
	drawnow
	
	%Create color map
	normalize = kMax-kMin;
	if(~normalize)
		normalize = 1;
	end
	cRs = abs((kMax - kVals)./normalize);
	clrs = [1-cRs,zeros(numKVals,1),cRs];
	
	%These are static variables, easier to just put them in a cell array
	plotPack = [par.domain,{par.plotVects},{par.useScatter},{figs},{clrs},{par.livePlot}];
	
	%*Data collection*
	[omg,eigVects,r] = par.getDataFunc(numKVals,kVals,logisParams,numPoints,par.T,par.dT,plotPack);
	
	[oM,oMInd] = max(real(omg));
	disp(['Max -omega^2 found at k = ' num2str(kVals(oMInd))])
	disp(['Largest value for -omega^2: ' num2str(oM)])
		
	PlotResults(r,eigVects,omg,kVals,plotPack{1:end-1})
	
	figstoprint = {figs{1:3},figs{6}};
	fignames = {'rho1','s1','phi1','omg'};
	foldername = sprintf('Output/numKvals%dkMin%0.3fkMax%0.3f/numPoints%d/',numKVals,kMin,kMax,numPoints);
	mkdirifnotexist(foldername);
	
	for i=1:numel(figstoprint)
		set(groot,'currentfigure',figstoprint{i});
		set(gcf,'color','w');
		export_fig(sprintf('%s%s.png', foldername, fignames{i}), '-r600', '-png');
	end
	
end

