function BilogisDemo
	
	small = 0.0001;
	c3 = small;
	r = (0:0.001:1);
	sz1 = numel(r);
	y1 = bilogis(r,c3,0,1,0,1);
	cr = [0,0,1];
	f1=figure('Position', [100, 200, 400, 300]);
	p1 = plot(r,y1,'color',cr);
	axis([0,1,0,1])
	title('Consistent curves between any two real subsets')
	
	y2 = halfbilogis(r,c3,0,1,0,1);
	f2=figure('Position', [600, 200, 400, 300]);
	p2 = plot(r,y2,'color',cr);
	axis([0,1,0,1])
	
	drawnow
	numiters = 100000;
	maxc3 = 30;
	
	maxlo1 = 10;
	maxhi1 = 30;
	maxnewlo1 = 3000;
	maxnewhi1 = 5000;
	
	for i=1:numiters
		c3 = maxc3*sin(2*pi/500*i)+small;
		
		lo1 = maxlo1*cos(2*pi/600*i);
		hi1 = maxhi1*sin(2*pi/200*i) + 41;
		
		newlo1 = maxnewlo1*sin(2*pi/100*i);
		newhi1 = maxnewhi1*cos(2*pi/300*i) + 9001;
		
		crr = abs(c3/(maxc3+small));
		cr = [crr,0,1-crr];
		r1 = linspace(lo1,hi1,sz1);
		if(c3>0)
			
			y1 = bilogis(r1,c3,lo1,hi1,newlo1,newhi1);
			y2 = halfbilogis(r,c3,0,1,0,1);
		else
			y1 = bilogit(r1,c3,lo1,hi1,newlo1,newhi1);
			y2 = halfbilogit(r,c3,0,1,0,1);
		end
		set(groot,'currentfigure',f1)
		axis([lo1,hi1,newlo1,newhi1])
		set(p1,'xdata',r1);
		set(p1,'ydata',y1);
		set(p1,'color',cr);
		set(p2,'ydata',y2);
		set(p2,'color',cr);
		drawnow
	end
	
	
end

