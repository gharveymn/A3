function Optimize
	x = [92.118,170.69,9.2188e+08,1000];
	delx2 = 1000;
	numpts = 1000;
	difs = ones(numpts,1);
	xs = zeros(numpts,4);
	for i = 1:numpts
		disp(i)
		difs(i) = OptFunc(x);
		xs(i,:) = x;
		
		x(1) = 100.*rand + 10;
		x(2) = delx2.*rand + 1;
		x(3) = 1e9*rand + 1e6;
	end
	M = horzcat(xs,difs);
	%N = csvread('Output/optData.dat');
	%M = cat(1,N,M(2:end,:));
	dlmwrite('Output/optData.txt',M,'precision',16)
end

