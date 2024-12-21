function mu0 = getMu0(a, r, s, G)

u0 = a(1:r);
b0 = a(r+1:r+s);
c0 = a(r+s+1:end);

mu0 = [];
flags = [];
for i = 1:s
    mu_init = 1;
    C1 = b0(i) - G(i);
    C2 = C1^3 - c0(i);
    options = optimoptions('fsolve', 'Display', 'none');
    [mu0_i,~,EXITFLAG] = fsolve(@(mu) equation_to_solve(mu, C1, C2), mu_init, options);
    flags = [flags EXITFLAG];
    mu0 = [mu0; mu0_i];
end

% % check
% for i =1:s
%     C1 = b0(i) - G(i);
%     C2 = C1^3 - c0(i);
%     equation_to_solve(mu0(i), C1, C2)
% end

function value = equation_to_solve(mu, C1, C2)

value = mu^3 - (abs(C1 - mu))^3 + C2;