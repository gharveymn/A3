function v=logis(v,x1,x2)
	m = x2-x1;
	sq = 10/m;
	v = 2*m.*(1./(1+exp(-sq*(v-x2)))) + x1;
end
