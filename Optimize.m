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

