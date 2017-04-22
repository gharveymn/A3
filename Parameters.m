function par=Parameters
	
	%Data analysis parameters
	par.plotVects = true;
	par.useScatter = true;
	par.livePlot = true;
	
	
	%Mathematical constraints
	par.curveParam = 500.597253981413900;%Higher implies higher slope
	par.rMin = 1e-4;
	par.rMax = 1e6;
	par.domain = {par.rMin,par.rMax};
	par.transformeddomain = {0,1.152811766819949e12};
	par.T = @(x1,x2,x3,x4,x5,x6) exp(x1);		%T for Transformation
	par.dT = @(xd,x1,x2,x3,x4,x5,x6) (-1).^(xd).*2.^(xd-1)./x1.^(xd+1);
	%par.T = @halfbilogit;		%T for Transformation
	%par.dT = @dhalfbilogit;
	
	
	%Data collection parameters
	par.numPoints = 1000;
	par.kMin = 0.001;
	par.kMax = 4;
	par.numKVals = 100;
	par.numEigs = 1;
	
	
end