function par=Parameters
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = true;
	par.livePlot = true;
	
	
	%Mathematical constraints
	par.curveParam = 50;%Higher implies higher slope
	par.rMin = -1000;
	par.rMax = 1000;
	par.domain = {par.rMin,par.rMax};
	par.transformeddomain = {-1e12,1e12};
	%par.T = @(x1,x2,x3,x4,x5,x6) exp(x1);		%T for Transformation
	%par.dT = @(xd,x1,x2,x3,x4,x5,x6) (-1).^(xd).*2.^(xd-1)./x1.^(xd+1);
	par.T = @bilogit;		%T for Transformation
	par.dT = @dbilogis;
	
	
	%Data collection parameters
	par.numPoints = 1000;
	par.kMin = 0.1;
	par.kMax = 0.5;
	par.numKVals = 10;
	par.numEigs = 1;
	
	
end