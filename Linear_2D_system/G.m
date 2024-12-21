function value = G(lambda,u,x1,N,m,r,Q,c,Q_c,control_bound)

% state constraints
value = [];
for k = 2:N+1 % recall that at k=1 x1 is specified
    xk = sol_diff_eqn(k, x1, u, m, r);
    for j = 1:size(c,2) % for each ellipsoid
        value = [value; -(xk - c(:,j))'*Q(:,2*(j-1)+1:2*(j-1)+2)*(xk - c(:,j)) + lambda];
    end
end

% control constraint:
if(~isempty(Q_c))
    u_parsed = parse_u(u, m, r);
    for k = 1:r/m
        value = [value; u_parsed(:,k)'*Q_c*u_parsed(:,k) - control_bound^2];
    end
end

