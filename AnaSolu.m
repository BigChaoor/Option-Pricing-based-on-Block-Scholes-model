function G = AnaSolu(mu1,mul2,sigma0, tao)
r = 0; y = 0.019;
q1 = (-y - 0.5*mu1*tao.^2 -1/3*mul2*tao.^3 - (0.5*sigma0^2+r)*tao)/sqrt(2*sigma0^2*tao);
q2 = (-y - 0.5*mu1*tao.^2 -1/3*mul2*tao.^3+ (0.5*sigma0^2-r)*tao)/sqrt(2*sigma0^2*tao);
G = 0.5*exp(-r*tao).*(exp(y+0.5*mu1*tao.^2+1/3*mul2*tao.^3+r*tao).*erfc(q1) - erfc(q2));

