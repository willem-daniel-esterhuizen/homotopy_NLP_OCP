function [direction_to_traverse_curve, t] = getTangent_QR(A, direction_to_traverse_curve)

% See Algower p.29
[Q,R] = qr(A');
t = Q(:,end);

if(strcmp(direction_to_traverse_curve, 'unset')) % i.e., the direction along which the curve must be tracked it not yet known
    if(t(1) > 0)
        if(det(Q)*det(R(1:end-1,:)) > 0) 
            direction_to_traverse_curve = 'positive';
        else
            direction_to_traverse_curve = 'negative';
        end
    else
        if(det(Q)*det(R(1:end-1,:)) > 0) 
            direction_to_traverse_curve = 'negative';
        else
            direction_to_traverse_curve = 'positive';
        end
    end
end

% This sign flip sometime does not work if det(R(1:end-1,:)) = inf
if(strcmp(direction_to_traverse_curve, 'positive'))
    if(det(Q)*det(R(1:end-1,:)) < 0)
        t = -t;
    end
else
    if(det(Q)*det(R(1:end-1,:)) > 0)
        t = -t;
    end
end


