function par=getPar
	par.plotImg = true;
	par.MP = false;
	par.useScatter = true;
	par.curveParam = 15;
	par.domain = {1e-4,200};
	par.transformeddomain = {1e1,1e9};
	par.numPoints = 1000;
	par.kMin = 0.2;
	par.kMax = 0.5;
	par.numKVals = 100;
end