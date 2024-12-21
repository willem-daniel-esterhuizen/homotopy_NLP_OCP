function plot_ellipsoids(Q, c)

for j = 1:size(c,2) % for each ellipsoid
    plot_ellipsoid_2D(Q(:,2*(j-1)+1:2*(j-1)+2), c(:,j));
end