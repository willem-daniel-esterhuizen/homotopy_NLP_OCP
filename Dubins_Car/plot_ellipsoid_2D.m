function plot_ellipsoid_2D(Q, centre)


% try : A_sqrt_inv * unit_sphere + centre

% The eigenvalues of positive definite Q satisfy x = 1/a^2, y = 1/b^2
% a, b semi-axes
E = eig(Q);

xc = centre(1);
yc = centre(2);

a = 1/sqrt(E(1));
b = 1/sqrt(E(2));

 % range to plot over
%------------------------------------
N = 50;
theta = 0:1/N:2*pi+1/N;

% Parametric equation of the ellipse
%----------------------------------------
state(1,:) = a*cos(theta); 
state(2,:) = b*sin(theta);

% Coordinate transform (since your ellipse is axis aligned)
%----------------------------------------
X = state;
X(1,:) = X(1,:) + xc;
X(2,:) = X(2,:) + yc;

% Plot
%----------------------------------------
plot(X(1,:),X(2,:),'k');
hold on;
axis equal;
grid;