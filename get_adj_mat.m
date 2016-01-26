function adj_mat = get_adj_mat( E, n )

% transform edges to adjacent matrix
adj_mat = zeros(n);
for i = 1:size(E, 1)
    adj_mat(E(i, 1), E(i, 2)) = 1;
    adj_mat(E(i, 2), E(i, 1)) = 1;
end

end