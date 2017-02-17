function [init_sets, labeled_cells] = SRG_graph( init_sets, cell_log_intensity, cell_area, n, adj_mat, invalid )
% graph-based seeded region growing
%
% Input variables:
%
% init_sets: (cell) the initial growing region sets represented by indices
% of voronoi cells
% cell_log_intensity: the log intensity of voronoi cells (corresponding to 
% photons)
% cell_area: the area of voronoi cells
% n: the number of voronoi cells
% adj_mat: the adjacent matrix that shows the connectivity of voronoi cells
% invalid: indices of invalid voronoi cells; invalid is due to that the
% area is zero/the voronoi cell is on the boundary
%
% Output variables:
%
% init_sets: (cell) the region sets after finish growing
% labeled_cells: the voronoi cells that have been assigned to one of the
% region sets

% initialize growing region sets
m = length(init_sets);
region_log_intensity = zeros(1, m);
region_area = zeros(1, m);
labeled_cells = [];
% neighboring voronoi cells of each growing region set
nbr_cells = cell(m, 1);
for j = 1:m
    % since init_sets{j} is a row vector
    % update labeled_cells
    labeled_cells = [labeled_cells init_sets{j}];
    region_area(j) = sum(cell_area(init_sets{j}));
    region_log_intensity(j) = log(sum(cell_area(init_sets{j}).*...
        exp(cell_log_intensity(init_sets{j})))/region_area(j));
end

% compute sets of neighbors for each growing region set
for j = 1:m
    nbr_cells{j} = cmp_nbr( init_sets{j}, adj_mat, n, [labeled_cells invalid] );
end

% get the number of voronoi cells that have been assigned to the initial 
% growing region sets
n_init = length(labeled_cells);
% get the number of invalid voronoi cells
n_invalid = length(invalid);
n_remain = n-n_init-n_invalid;

% start the iteration to assign the remaining voronoi cells 
for i = 1:n_remain
    
    if mod(i, 100)==0
        i
    end
    
    % select the pair of a growing region set and a neighboring voronoi 
    % cell with the smallest value of the criterion
    min_c = realmax;
    % j refers to regions, and k refers to cells
    for j = 1:m
        for k = nbr_cells{j}
            % the criterion is the abs difference between the log 
            % intensities
            c = abs(cell_log_intensity(k)-region_log_intensity(j));
            if c<min_c
                min_c = c;
                min_j = j;
                min_k = k;
            end
        end
    end
    
    % update the selected growing region set
    labeled_cells = [labeled_cells min_k];
    init_sets{min_j} = [init_sets{min_j} min_k];
    region_log_intensity(min_j) = log((exp(region_log_intensity(min_j))*...
        region_area(min_j)+exp(cell_log_intensity(min_k))*cell_area(min_k))/...
        (region_area(min_j)+cell_area(min_k)));
    region_area(min_j) = region_area(min_j)+cell_area(min_k);
    % add the set of neighbors for cell min_k to that of region min_j
    nbr_cells{min_j} = [nbr_cells{min_j} cmp_nbr( min_k, adj_mat, n, [labeled_cells invalid] )];
    % delete cell min_k from the set of neighbors for each growing region set
    for j = 1:m
        nbr_cells{j} = setdiff(nbr_cells{j}, min_k);
    end
    
end
        
end