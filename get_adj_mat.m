function adj_mat = get_adj_mat( E, n )
% transform edges to adjacent matrix
% E is obtained from function edges
% adj_mat: adjacent matrix, 1 represents that there is an edge

adj_mat = false(n);
for i = 1:size(E, 1)
    adj_mat(E(i, 1), E(i, 2)) = true;
    adj_mat(E(i, 2), E(i, 1)) = true;
end

end