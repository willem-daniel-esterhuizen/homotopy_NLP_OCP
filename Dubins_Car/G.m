function value = G(lambda,u,x1,N,m,r,Q,c,Q_c,min_car_turning_radius,v)

%state constraints
value = [];
for k = 2:N+1 % recall that at k=1 x1 is specified
    xk = sol_diff_eqn(k, x1, u, v);
    for j = 1:size(c,2) % for each ellipsoid
        value = [value; -(xk(1:2) - c(:,j))'*Q(:,2*(j-1)+1:2*(j-1)+2)*(xk(1:2) - c(:,j)) + lambda];
    end
end

% control constraint:
if(~isempty(Q_c))
    for k = 1:r/m
        value = [value; u'*Q_c*u - (1/min_car_turning_radius)^2];
    end
end

