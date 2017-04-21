function par=Parameters
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = true;
	
	
	%Mathematical constraints
	par.curveParam = 1000;%Higher implies higher slope
	par.domain = {0.001,1000};
	par.transformeddomain = {1e-16,1e25};
	par.T = @halfbilogit;		%T for Transformation
	par.dT = @dhalfbilogit;
	
	
	%Data collection parameters
	par.numPoints = 10000;
	par.kMin = 0.01;
	par.kMax = 4;
	par.numKVals = 100;
	par.numEigs = 1;
	
	
end