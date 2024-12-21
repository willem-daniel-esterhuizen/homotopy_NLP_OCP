function rho = Rho(lambda, u, mu, a, N, r, s, x1, Jacobian_u_J, Jacobian_u_G, m, Q, c, Q_c,min_car_turning_radius,v)

u0 = a(1:r);
b0 = a(r+1:r+s);
c0 = a(r+s+1:end);

G_sym = G(lambda,u,x1,N,m,r,Q,c,Q_c, min_car_turning_radius,v);

rho = lambda*(Jacobian_u_J' + Jacobian_u_G'*mu) + (1 - lambda)*(u - u0);
for i = 1:s
    rho = [rho; 
           mu(i)^3 - (abs((1 - lambda)*b0(i) - G_sym(i) - mu(i)))^3 + ( (1 - lambda)*b0(i) - G_sym(i) )^3 - (1 - lambda)*c0(i)];
end