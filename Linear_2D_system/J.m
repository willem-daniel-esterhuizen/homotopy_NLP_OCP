function value = J(u, N, x1, x_target, m, r)
% the cost function

% sum_{k=1}^{N} 0.5 norm(u_k)^2 + 0.5*norm(x_{N + 1} - x_target)^2
u_parsed = parse_u(u, m, r);
value = 0;
for k=1:N
    value = value + 0.5*norm(u_parsed(:,k))^2;
end

x_N_plus_1 = sol_diff_eqn(N+1, x1, u, m, r);
value = value + 0.5*norm(x_N_plus_1 - x_target)^2;