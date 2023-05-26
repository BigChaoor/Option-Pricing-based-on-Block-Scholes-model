tao = 1:0.05:2;
n = length(tao);
SigmaEps = 1E-2; 
SigmaTheta = [0.01, 0; 0, 0.01];  
G = AnalyticalSolu(1,1,tao);
T = G.*(1 + SigmaEps*randn(1,n));
%%%%%%%%%-------%%%%%%%%%
y = -15:30/199:15;
F = Crank_Nicholson(1,1,y);
F = F(91:110);
m = length(F);
Y = F.*(1 + SigmaEps*randn(m,1));

%% Metropolis-Hastings Algorithm
theta1 = [0.5,0.5];
theta2 = [4,4];
k_star = 5000;
sample1 = zeros(k_star,2);  
sample2 = zeros(k_star,2);
for k = 1:5000
    theta_Skim = mvnrnd(theta2,SigmaTheta);
    G_Skim = AnalyticalSolu(theta_Skim(1),theta_Skim(2),tao);
    Gk = AnalyticalSolu(theta2(1),theta2(2),tao);
    b0 = (T - Gk)*(T - Gk)'-(T - G_Skim)*(T - G_Skim)';
    alpha0 = min(exp(b0/(2*SigmaEps^2)),1);
    %%%%%%%%------------------%%%%%%%%
    F_Skim = Crank_Nicholson(theta_Skim(1),theta_Skim(2),y);
    Fk = Crank_Nicholson(theta1(1),theta1(2),y);
    b1 = (Y - Fk(91:110))'*(Y - Fk(91:110)) - (Y - F_Skim(91:110))'*(Y - F_Skim(91:110));
    alpha1 = min(exp(b1/(2*SigmaEps^2)),1);  
    U = rand;
    if U<alpha0
        theta2 = theta_Skim;
    end
    if U<alpha1
        theta1 = theta_Skim;
    end
    sample1(k,:) = theta1;
    sample2(k,:) = theta2;
end
fprintf('theta0=%2.15f\n',theta2)
fprintf('theta1=%2.15f\n',theta1)

%% plot the figures of Time independent case 
x = 1:k_star;
ture = ones(1,k_star);
figure(1)
plot(x,sample1(:,1),'b','linewidth',1)
hold on
plot(x,true,'r','linewidth',1)
legend('MCMC sapmling of the parameter mu1','true value of mu1');
ylim([min(sample1(:,1))-0.2,max(sample1(:,1))+0.2])
title('Time independent case')

figure(2)
plot(x,sample1(:,2),'b','linewidth',1)
hold on
plot(x,ture,'r','linewidth',1)
legend('MCMC sapmling of the parameter sigma0','true value of sigma0');
ylim([min(sample1(:,2))-0.2,max(sample1(:,2))+0.2])
title('Time independent case')
%%%%%%%%%%--------%%%%%%%%

%% plot the figures of Time dependent case 
figure(3)
plot(x,sample2(:,1),'b','linewidth',1)
hold on
plot(x,ture,'r','linewidth',1)
legend('MCMC sapmling of the parameter mu1','true value of mu1');
ylim([min(sample2(:,1))-0.2,max(sample2(:,1))+0.2])
title('Time dependent case')

figure(4)
plot(x,sample2(:,2),'b','linewidth',1)
hold on
plot(x,ture,'r','linewidth',1)
legend('MCMC sapmling of the parameter sigma0','true value of sigma0');
ylim([min(sample2(:,2))-0.2,max(sample2(:,2))+0.2])
title('Time dependent case')

%% L-M Algorithm
% theta_LM1 = lsqcurvefit(@CNSolu,theta0,y(90:110)',Y)
theta_LM2 = lsqcurvefit(@ExactSolu,[4,4],tao,T)