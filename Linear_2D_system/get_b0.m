function b0_min = get_b0(u0,x1,N,m,r,s,Q,c,Q_c,control_bound)

value = G(0,u0,x1,N,m,r,Q,c,Q_c,control_bound);

b0_init = ones(s,1);
options = optimoptions('fmincon', 'Display', 'none');
b0_min = fmincon(@(b0) ones(1,s)*b0, b0_init, [],[],[],[],[],[], @(b0) constraints(b0, value), options);

% if b0 is too small, then (b0_i - G_i(0,u0))^3 might be very small, which
% then buggers up the computation of mu0...
b0_min = max(ones(s,1), b0_min); % recall b0 > 0.

function [c,ceq] = constraints(b0, value)

ceq = [];
c = value - b0;

