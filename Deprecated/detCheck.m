sz = 20;
maxVal = 10;


A = diag(maxVal.*rand(1,sz));
B = diag(maxVal.*rand(1,sz));
C = diag(maxVal.*rand(1,sz));
D = diag(maxVal.*rand(1,sz));
E = diag(maxVal.*rand(1,sz));
F = diag(maxVal.*rand(1,sz));
G = diag(maxVal.*rand(1,sz));
H = diag(maxVal.*rand(1,sz));
K = diag(maxVal.*rand(1,sz));

T = [A B C
	D E F
	G H K];

P =  [A B C
	D E F
	G zeros(sz) K];

Q = [eye(sz) zeros(sz) zeros(sz)
	zeros(sz) zeros(sz) zeros(sz)
	zeros(sz) zeros(sz) zeros(sz)];

dettest = 1;

for n = 1:size(A,1)
	dettest = dettest*1/A(n,n)*((A(n,n)*K(n,n) - C(n,n)*G(n,n))*(A(n,n)*E(n,n) - B(n,n)*D(n,n)) ...
	- (A(n,n)*F(n,n) - C(n,n)*D(n,n))*(A(n,n)*H(n,n) - B(n,n)*G(n,n)));
end


realDet = det(T);

disp(['myDet: ' num2str(dettest)]);
disp(['realDet: ' num2str(realDet)]);

eigens = zeros(0,1);
for n = 1:size(A,1)
	eigens = cat(1,eigens,roots([1,-(2*A(n,n) - B(n,n)*(D(n,n) - 1) - C(n,n)),(A(n,n)-C(n,n))*(A(n,n)-B(n,n)*D(n,n))]));
end

scatter(real(eigens),imag(eigens),[],'r','.');
hold on
opts = struct('p', sz);
eins = eigs(P,Q,18,'sm',opts);
scatter(real(eins),imag(eins),[],'b','.');

