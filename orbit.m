function orbit
	a = 4;
	b = 1;
	eps = sqrt(1 - (b/a)^2);
	
	l = 1;
	gamma = 1;
	mu = 1;
	
	c = l^2/(gamma*mu);
	
	phi = (0:0.01:10)';
	
	r = c./(1 + eps.*cos(phi));
	
	x = r.*cos(phi);
	y = r.*sin(phi);
	
	p1 = plot(x,y);
	axis([-35,5,-4,4])
	
	for i=0:size(x,1)
		set(p1,'XData',x(1:i));
		set(p1,'YData',y(1:i));
		drawnow;
	end
	
end

