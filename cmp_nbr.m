function nbr_cells = cmp_nbr( seed_set, adj_mat, n, labeled_cells )
% compute neighbors of a seed set

% exclude seed_set itself and labeled_cells since they cannot be neighbors
diff_set = setdiff(1:n, [seed_set labeled_cells]);

% compute sum of each row
nbr_cells = diff_set(sum(adj_mat(diff_set, seed_set), 2)>0);

end