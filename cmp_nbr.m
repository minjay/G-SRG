function nbr_set = cmp_nbr( seed_set, adj_mat, n, labeled_set )

% compute neighbors of a seed set
diff_set = setdiff(1:n, [seed_set labeled_set]);
nbr_set = [];
for i = diff_set
    for j = seed_set
        if adj_mat(i, j)==1
            nbr_set = [nbr_set i];
            break
        end
    end
end

end