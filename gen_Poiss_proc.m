function [flux, V, R] = gen_Poiss_proc(mu, V, R, max_t)

n = length(R);
n_p = poissrnd(mu);
cx = [];
cy = [];

for i = 1:n
    if isnan(n_p(i))
        continue
    end
    max_x = max(V(R{i}, 1));
    min_x = min(V(R{i}, 1));
    max_y = max(V(R{i}, 2));
    min_y = min(V(R{i}, 2));
    rx = rand(1, max_t).*(max_x-min_x)+min_x;
    ry = rand(1, max_t).*(max_y-min_y)+min_y;
    flag = inpolygon(rx, ry, V(R{i}, 1), V(R{i}, 2));
    tx = rx(flag==1);
    ty = ry(flag==1);
    cx = [cx tx(1:n_p(i))];
    cy = [cy ty(1:n_p(i))];
end

DT = delaunayTriangulation([cx' cy']);
[V, R] = voronoiDiagram(DT);
n = length(cx);
flux = zeros(n, 1);
area_cell = zeros(n, 1);
for i = 1:n
    area_cell(i) = polyarea(V(R{i}, 1), V(R{i}, 2));
    flux(i) = 1/area_cell(i);
end

% triplot(DT)
% hold on
% scatter(cx', cy', [], log(flux), 'o', 'filled')
% colorbar
% axis image
    