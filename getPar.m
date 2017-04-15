function par=getPar
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = true;
	
	
	%Mathematical constraints
	par.curveParam = 30;%Higher implies higher slope
	par.domain = {-1000,1000};
	par.transformeddomain = {-1e21,1e21};
	par.T = @bilogis;		%T for Transformation
	par.dT = @dbilogis;
	
	
	%Data collection parameters
	par.numPoints = 10000;
	par.kMin = 0.01;
	par.kMax = 0.5;
	par.numKVals = 5;
	par.numEigs = 1;
	
	
end