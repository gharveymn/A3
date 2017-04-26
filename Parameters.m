function par=Parameters
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = true;
	par.livePlot = true;
	
	
	%Mathematical constraints
	par.curveParam = 10;%Higher implies higher slope
	par.rMin = -10;
	par.rMax = 10;
	par.domain = {par.rMin,par.rMax};
	par.transformeddomain = {-1e6,1e6};		%For bilogis/bilogit
	par.T = @(x1,x2,x3,x4,x5,x6) exp(x1);		%T for Transformation
	par.dT = @(xd,x1,x2,x3,x4,x5,x6) (-1).^(xd+1)./x1.^(xd);
	%par.T = @halfbilogit;		%T for Transformation
	%par.dT = @dhalfbilogis;
	
	
	%Data collection parameters
	par.numPoints = 100;
	par.kMin = 0.1;
	par.kMax = 0.5;
	par.numKVals = 20;
	par.numEigs = 1;
	
	
end