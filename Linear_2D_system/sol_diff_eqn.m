function x_sol = sol_diff_eqn(M, x1, u, m, r)

% returns the state at k = M, starting from x1 at k = 1, with input
% sequence u(k), k = 1:N-1

u_parsed = parse_u(u, m, r);

x_sol = x1;
for k = 1:M-1
    x_sol = x_sol + u_parsed(:,k); % specify dynamics here % x_k+1 = x_k + u_k
end