function x_sol = sol_diff_eqn(M, x1, u, v)

% returns the state at k = M, starting from x1 at k = 1, with input
% sequence u(k), k = 1:M-1

h = 0.1;
x_sol = x1;
for k = 1:M-1
    x_sol = [x_sol(1) + v*h*cos(x_sol(3));
             x_sol(2) + v*h*sin(x_sol(3));
             x_sol(3) + h*u(k);];            
end