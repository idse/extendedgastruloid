function map = cmap_rand(n)

map(2:n,1:3) = hsv(double(n-1));
[~, order] = sort(rand(n-1,1));
map(2:n,1:3) = map(order+1,:);
