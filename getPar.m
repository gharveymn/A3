function par=getPar
	par.plotVects = true;
	par.useScatter = true;
	par.curveParam = 30;%Higher implies higher slope
	par.domain = {1e-16,1000};
	par.transformeddomain = {1e-16,1e21};
	par.numPoints = 10000;
	par.kMin = 0.01;
	par.kMax = 0.5;
	par.numKVals = 5;
	par.numEigs = 1;
end