function F = Crank_Nicholson(mu1,sigma0,y)
r = 0; N = length(y)-1;
DetaTao = 0.01;
DetaY = (y(end)-y(1))/N; 
tao_star = 0.55;
M = tao_star/DetaTao;
% y = -15:DetaY:15; 
a = -DetaTao/DetaY^2/4*(sigma0+DetaY*(-0.5*sigma0+mu1*y(2:N-1)));
b = DetaTao/DetaY^2/2*ones(1,N-1);
c = -DetaTao/DetaY^2/4*(sigma0-DetaY*(-0.5*sigma0+mu1*y(3:N)));
b_tilde = b + DetaTao*r;
A = diag(1+b) + diag(c,-1) + diag(a,1);
B = diag(1-b_tilde) + diag(-c,-1) + diag(-a,1);
U0 = max(exp(y(2:N))-1,0)';
% U = zeros(size(U0));
for i = 1:M
    U = A\(B*U0);
    U0 = U;
end
F = U;