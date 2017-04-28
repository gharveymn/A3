function Optimize
	x = [20,100,1e6,100];
	delx2 = 100;
	numpts = 20;
	difs = ones(numpts,1);
	xs = zeros(numpts,4);
	for i = 1:numpts
		disp(i)
		difs(i) = OptFunc(x);
		xs(i,:) = x;
		
		x(1) = 30.*rand + 5;
		x(2) = delx2.*rand + 1;
		x(3) = 1e9*rand + 1e3;
	end
	M = horzcat(xs,difs);
	%N = csvread('Output/optData.dat');
	%M = cat(1,N,M(2:end,:));
	dlmwrite('Output/optData.txt',M,'precision',16)
end

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

