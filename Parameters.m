function par=Parameters
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = false;
	par.livePlot = true;
	
	
	%Mathematical constraints
	par.numPoints = 100;
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
	par.numKVals = 1000;
	par.kMin = 0.209;
	par.kMax = 0.225;
	par.numEigs = 1;
	par.getDataFunc = @GetData;
	
	
end