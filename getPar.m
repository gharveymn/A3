function par=getPar
	par.plotImg = true;
	par.MP = false;
	par.useScatter = true;
	par.curveParam = 15;
	par.domain = {1e-4,10};
	par.transformeddomain = {1,1e12};
	par.numPoints = 100;
	par.kMin = 0.01;
	par.kMax = 0.4;
	par.numKVals = 10;
	par.figs = {figure(1),figure(2),figure(3),figure(4),figure(5),figure(6)};
end