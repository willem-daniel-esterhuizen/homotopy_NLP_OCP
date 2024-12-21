function value = J(u, N, x1, x_target, v)

%cost function

value = 0;
% for k=1:N
%     value = value + 0*norm(u(k))^2;
% end
% for k=2:N
%     xk = sol_diff_eqn(k, x1, u);
%     value = value + norm(xk(1:2) - x_target(1:2))^2;
% end

% % Minimise distance to target at the end:
% x_N_plus_1 = sol_diff_eqn(N+1, x1, u);
% value = value + norm(x_N_plus_1(1:2) - x_target(1:2))^2; % only care about position

% minimise arc length

% for k=2:N+1
%     x_k_minus_1 = sol_diff_eqn(k-1, x1, u);
%     x_k = sol_diff_eqn(k, x1, u);
%     value = value + norm(x_k_minus_1(1:2) - x_k(1:2))^2;
% end
x_N_plus_1 = sol_diff_eqn(N+1, x1, u, v);
value = value + norm(x_N_plus_1(1:2) - x_target(1:2))^2;