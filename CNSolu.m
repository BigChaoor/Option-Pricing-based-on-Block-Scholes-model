function F = CNSolu(d,y)
r = 0; N = length(y);
DetaTao = 0.01;
DetaY = (y(end)-y(1))/N; 
tao_star = 0.55;
M = tao_star/DetaTao;
% y = -15:DetaY:15; 
a = -DetaTao/DetaY^2/4*(d(2)+DetaY*(-0.5*d(2)+d(1)*y(1:N-1)));
b = DetaTao/DetaY^2/2*ones(1,N);
c = -DetaTao/DetaY^2/4*(d(2)-DetaY*(-0.5*d(2)+d(1)*y(2:N)));
b_tilde = b + DetaTao*r;
A = diag(1+b) + diag(c,-1) + diag(a,1); 
B = diag(1-b_tilde) + diag(-c,-1) + diag(-a,1);
U0 = max(exp(y)-1,0)';
U = zeros(size(U0));
for i = 1:M
    U = A\(B*U0);
    U0 = U;
end
F = U;