y = -15:30/199:15;
F = Crank_Nicholson(1,1,y);
F = F(91:110);
m = length(F);
SigmaEps = 1E-2; 
Y = F.*(1 + SigmaEps*randn(m,1));
%%%%%%%%%%%%%%
SigmaTheta = [0.01, 0; 0, 0.01];  
k_star = 5000;
sample1 = zeros(k_star,2);
theta1 = [0.5,0.5];
for k = 1:5000
    theta_Skim = mvnrnd(theta1,SigmaTheta);
    F_Skim = Crank_Nicholson(theta_Skim(1),theta_Skim(2),y);
    Fk = Crank_Nicholson(theta1(1),theta1(2),y);
    b1 = (Y - Fk(91:110))'*(Y - Fk(91:110)) - (Y - F_Skim(91:110))'*(Y - F_Skim(91:110));
    alpha1 = min(exp(b1/(2*SigmaEps^2)),1);  
    U = rand;
    if U<alpha1
        theta1 = theta_Skim;
    end
    sample1(k,:) = theta1;
end
x = 1:k_star;
figure(3)
plot(x,sample1(:,1),'b')
hold on
plot(x,ones(1,k_star),'r','linewidth',1)
legend('MCMC sapmling of the parameter mu1','true value of mu1');
% ylim([-0.1,1.25])
%%%%%%%
figure(4)
plot(x,sample1(:,2),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter sigma0','true value of sigma0');
ylim([-0.1,max(theta1(:,2))])