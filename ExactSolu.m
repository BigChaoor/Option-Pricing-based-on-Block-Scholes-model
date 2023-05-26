function z = ExactSolu(a, x)
r = 0; y = 0.019;
q1 = (-y - 0.5*a(1)*x.^2 - (0.5*a(2)^2+r)*x)/sqrt(2*a(2)^2*x);
q2 = (-y - 0.5*a(1)*x.^2 + (0.5*a(2)^2-r)*x)/sqrt(2*a(2)^2*x);
z = 0.5*exp(-r*x).*(exp(y+0.5*a(1)*x.^2+r*x).*erfc(q1) - erfc(q2));