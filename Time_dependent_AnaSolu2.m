%% parameter settings
tao = 1:0.05:2;
n = length(tao);
SigmaEps = 1E-2; 
G = AnalyticalSolu(1,1,tao);
T = G.*(1 + SigmaEps*randn(1,n));

%% MCMC algorithm via AnalyticalSolu2
SigmaTheta = [0.01, 0; 0, 0.01];  
k_star = 5000;
sample2 = zeros(k_star,2);
theta0 = [4,4];
for k = 1:k_star
    theta_Skim = mvnrnd(theta0,SigmaTheta);
    G_Skim = AnalyticalSolu(theta_Skim(1),theta_Skim(2),tao);
    Gk = AnalyticalSolu(theta0(1),theta0(2),tao);
    b0 = (T - Gk)*(T - Gk)'-(T - G_Skim)*(T - G_Skim)';
    alpha0 = min(exp(b0/(2*SigmaEps^2)),1);
    U = rand;
    if U<alpha0
        theta0 = theta_Skim;
    end
    sample2(k,:) = theta0;
end
%% figures
x = 1:k_star;
figure(1)
plot(x,sample2(:,1),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter mu1','true value of mu1');
% ylim([-0.1,1.25])
%%%%
figure(2)
plot(x,sample2(:,2),'b')
hold on
plot(x,ones(1,k_star),'r')
legend('MCMC sapmling of the parameter sigma0','true value of sigma0');
% ylim([-0.1,1.25])
%% L-M Algorithm
theta = lsqcurvefit(@ExactSolu,theta0,tao,T)
