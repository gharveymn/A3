function maxDif=OptFunc(xpar)
	
	%ARGS
	curveParam = xpar(1);
	%domain1    = -xpar(2);
	domain1	 = 0;
	domain2    = xpar(2);
	%tdomain1   = -xpar(3);
	tdomain1	 = 0;
	tdomain2   = xpar(3);
	numPoints  = xpar(4);
	
	%ASSIGNED
	%numPoints = 100;
	
	
	%Get parameters
	numKVals = 10;
	kMin = 0.20;
	kMax = 0.41;
	kVals = linspace(kMin,kMax,numKVals)';
	logisParams = [{curveParam},{domain1,domain2},{tdomain1,tdomain2}];
	sz = numPoints;
	
	plotPack = {0,0,false,0,0,0,false};
	
	%*Data collection*
	[omg,eigVects,r] = GetData(numKVals,kVals,logisParams,sz,@halfbilogit,@dhalfbilogis,plotPack);
	%abs(diff(omg))
	maxDif = max(abs(diff(omg)));
	if(~maxDif)
		maxDif = inf;
	end
	
	disp(xpar)
	%disp(omg)
	disp(maxDif)
	
end

