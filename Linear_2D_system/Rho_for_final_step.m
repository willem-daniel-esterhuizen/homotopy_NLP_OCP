function rho = Rho_for_final_step(y, a, r, s, x1, Jacobian_u_J, Jacobian_u_G, lambda, u, mu, m, Q, c)

lambda_val = 1;
u_val = y(1:r);
mu_val = y(r+1:end);

Jacobian_u_J = subs(Jacobian_u_J, [lambda; u; mu], [lambda_val; u_val; mu_val]);
Jacobian_u_G = subs(Jacobian_u_G, [lambda; u; mu], [lambda_val; u_val; mu_val]);
Jacobian_u_J = eval(Jacobian_u_J);
Jacobian_u_G = eval(Jacobian_u_G);

rho = Rho(lambda_val, u_val, mu_val, a, r, s, x1, Jacobian_u_J, Jacobian_u_G, m, Q, c);
