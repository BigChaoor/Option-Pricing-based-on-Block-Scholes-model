tao = 1:0.05:2;
n = length(tao);
SigmaEps = 0.001; 
G = AnalyticalSolu2(1,1,1,tao);
T = G.*(1 + 0.001*randn(1,n));
%%%%%%%%%%%%%%%%%%%%
SigmaTheta = diag([0.01,0.01,0.01]);  
W =  diag([-1,-1,-1]);
k_star = 1000000;
sample2 = zeros(k_star,3);
theta0 = [0.1,0.1,0.3];
for k = 1:k_star
    s = mvnrnd([0,0,0],SigmaTheta);
    theta_Skim =sqrt(1-SigmaEps^2).*theta0+SigmaEps*s;
    if theta_Skim(3)>0
    G_Skim = AnalyticalSolu2(theta_Skim(1),theta_Skim(2),theta_Skim(3),tao);
    Gk = AnalyticalSolu2(theta0(1),theta0(2),theta0(3),tao);
    b0 = (T - Gk)*(T - Gk)'-(T - G_Skim)*(T - G_Skim)'+0*SigmaEps^2.*(theta0*W*theta0'-theta_Skim*W*theta_Skim');
    alpha0 = min(exp(b0/(2*SigmaEps^2)),1);
    U = rand;
    if U<alpha0
        theta0 = theta_Skim;
    end
    end
    sample2(k,:) = theta0;
   
end
x = 1:k_star;
figure(1)
plot(x,sample2(:,1),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter mu1','true value of mu1');
ylim([-0.1,1.25])
%%%%
figure(2)
plot(x,sample2(:,2),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter mu12','true value of mu12');
ylim([-0.1,1.25])
%%%%

%%%%
figure(3)
plot(x,sample2(:,3),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter sigma0','true value of sigma0');
ylim([-0.1,1.25])