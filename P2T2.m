function P2T2
	
	rMax = 20;
	%numPts = 1000000;
	
	h = 0.00001;
	rinit = (0.00001:h:rMax)';
	
	%squeeze the domain
	r = logis(rinit,0,rMax);
	%r = rinit;
	rInv = 1./r;
	p0Curr = exp(-1./10.*r.^2);
	%p0Curr = 10000000.*sin(1000.*r).*cos(300*r).*tan(40.*r);
	%p0Curr(end) = 0;
	p0Last = zeros(size(p0Curr,1),1);
	phi0 = p0Last;
	
	%p0Curr = 1./(1 + r.^2./8).^2;
	%p0Last = zeros(size(p0Curr,1)+2,1);
	%phi0 = -log(p0Curr);
	
	eps = 10^-25;
	
	%p0Comp = 1./(1 + r.^2./8).^2;
	%phi0Comp = -log(p0Comp);
	
	f1 = figure(1);
	clf
	hold on
	p1=plot(r,p0Curr);
	p2=plot(r,phi0);
	%p3=plot(r,p0Comp);
	%p4=plot(r,phi0Comp);
	%legend([p1,p2,p3,p4],'p0Curr','phi0','p0Comp','phi0Comp');
	%axis([0,20,0,1]);
	drawnow;
	
	difVec = (p0Last - p0Curr).^2;
	
	f2 = figure(2);
	p3 = plot(rinit,difVec);
	title('Numerical Change')
	drawnow;
	
	i = 0;
	
	d = eps+1;
	
	while(d > eps)
		
		p0Last = p0Curr;
		
		%Internal method
		phi0 = cumtrapz(r,rInv.*cumtrapz(r,r.*p0Curr));
		
		%Vectorized Trapezoidal
		%b1 = r.*p0Curr;
		%b2 = rInv.*cumsum((b1 + [diff(b1);b1(end)-b1(end-1)]./2).*h);
		%phi0 = cumsum((b2 + [diff(b2);b2(end)-b2(end-1)]./2).*h);
		
		%Rectangle Method
		%phi0 = cumsum(rInv.*cumsum(r.*p0Curr.*h).*h);
		
		p0Curr = exp(-phi0);
		%rInv = rInv(1:end-2);
		%r = r(1:end-2);
		
		set(groot,'currentfigure',f1);
		%set(p1,'XData',r);
		set(p1,'YData',p0Curr);
		%set(p2,'XData',r);
		set(p2,'YData',phi0);
		drawnow;
		i = i +1;
		
		difVec = (p0Last - p0Curr).^2;
		set(groot,'currentfigure',f2);
		set(p3,'YData',difVec);
		drawnow;
		
		
		d = sum(difVec);
		%d = sum((p0Last(1:10:end) - p0Curr(1:10:end)).^2);
		disp(d)
	end
	disp(i)
	
	difVec = (p0Last - p0Curr).^2;
	
	figure(2)
	plot(r,difVec);
	
	
end

function v=logit(v,m)
	v = -log(2*m./(v+m)-1);
end

