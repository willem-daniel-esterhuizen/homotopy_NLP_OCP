function t = getTangent_QR(A)

% See Algower p.29

[Q,R] = qr(A');
t = Q(:,end);

if(det(Q)*det(R(1:end-1,:)) < 0)
    t = -t;
end

