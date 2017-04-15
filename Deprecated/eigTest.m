function eigTest
%syms x
maxVal = 1;
pts = zeros(0,2);

for i=1:10000
rng = maxVal.*rand(1,9);
A = rng(1);
B = rng(2);
C = rng(3);
D = rng(4);
E = rng(5);
F = rng(6);
G = rng(7);
H = rng(8);
K = rng(9);
%lam = linspace(-10,10);

alpha = A*K-C*G;
beta = A+K;
gamma = A*E-B*D;
delta = A+E;

epsilon = A*F-C*D;
eta = A*H-B*G;

%y=lam.^4-(beta+delta).*lam.^3+(alpha+beta*delta-F*H).*lam.^2-(alpha*delta+beta*gamma-epsilon*H-eta*F).*lam...
%												+alpha*gamma-epsilon*eta;
%{											
plot(lam,y)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
%}
%p1 = x.^4-(beta+delta).*x.^3+(alpha+beta*delta-F*H).*x.^2-(alpha*delta+beta*gamma-epsilon*H-eta*F).*x...
%												+alpha*gamma-epsilon*eta;
										
p = [1,-(beta+delta),(alpha+beta*delta-F*H),-(alpha*delta+beta*gamma-epsilon*H-eta*F),alpha*gamma-epsilon*eta];

S = roots(p);

re = real(S);
im = imag(S);

re = re(re(:) ~= A);
im = im(re(:) ~= A);

pts = cat(1,pts,[re im]);

end

scatter(pts(:,1),pts(:,2),[],'.');
end
