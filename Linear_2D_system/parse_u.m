function u_parsed = parse_u(u, m, r)

% takes the r x 1 vector and returns an m x N vector
% r = N x m
u_parsed = [];
for i = 1:r/m
    u_parsed = [u_parsed, [u(i*m - 1);u(i*m)]];
end