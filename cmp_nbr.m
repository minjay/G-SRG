function nbr_cells = cmp_nbr( seed_set, adj_mat, n, labeled_cells )
% compute neighbors of a seed set

% exclude seed_set itself and labeled_cells since they cannot be neighbors
diff_set = setdiff(1:n, [seed_set labeled_cells]);

nbr_cells = [];
for i = diff_set
    for j = seed_set
        if adj_mat(i, j)==1
            nbr_cells = [nbr_cells i];
            break
        end
    end
end

end