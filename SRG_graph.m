function [init_sets, labeled_cells] = SRG_graph( init_sets, cell_log_intensity, cell_area, n, adj_mat, invalid, verbose, print_n )
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
% verbose: whether print out the status
% print_n: print every n iterations
%
% Output variables:
%
% init_sets: (cell) the region sets after finish growing
% labeled_cells: the voronoi cells that have been assigned to one of the
% region sets

if nargin==6
    verbose = false;
end

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

% min_c_pairs contains the min diff and the corresponding cell index
min_c_pairs = zeros(m, 2);
for j = 1:m
    % if region j has neighbors
    if ~isempty(nbr_cells{j})
        [min_c, min_index] = min(abs(cell_log_intensity(nbr_cells{j})-region_log_intensity(j)));
        min_c_pairs(j, 1) = min_c;
        min_c_pairs(j, 2) = nbr_cells{j}(min_index);
    else
        min_c_pairs(j, 1) = realmax;
        min_c_pairs(j, 2) = 0;
    end
end

% get the number of voronoi cells that have been assigned to the initial 
% growing region sets
n_init = length(labeled_cells);
% get the number of invalid voronoi cells
n_invalid = length(invalid);
n_remain = n-n_init-n_invalid;

% start the iteration to assign the remaining voronoi cells 
for i = 1:n_remain
    
    if verbose
        if mod(i, print_n)==0
            disp('Iteration ' + str(i) + '...')
        end
    end
    
    % select the pair of a growing region set and a neighboring voronoi 
    % cell with the smallest value of the criterion
    % j refers to regions, and k refers to cells
    [~, min_j] = min(min_c_pairs(:, 1));
    min_k = min_c_pairs(min_j, 2);
    
    % update the selected growing region set
    labeled_cells = [labeled_cells min_k];
    init_sets{min_j} = [init_sets{min_j} min_k];
    region_log_intensity(min_j) = log((exp(region_log_intensity(min_j))*...
        region_area(min_j)+exp(cell_log_intensity(min_k))*cell_area(min_k))/...
        (region_area(min_j)+cell_area(min_k)));
    region_area(min_j) = region_area(min_j)+cell_area(min_k);
    % add the set of neighbors for cell min_k to that of region min_j
    % make sure the uniqueness of the elements
    nbr_cells{min_j} = unique([nbr_cells{min_j} cmp_nbr( min_k, adj_mat, n, [labeled_cells invalid] )]);
   
    for j = 1:m
        if ismember(min_k, nbr_cells{j})
            % delete cell min_k from the set of neighbors for each growing region set
            nbr_cells{j} = setdiff(nbr_cells{j}, min_k);
            if min_c_pairs(j, 2)==min_k
                % update min_c_pairs for region j
                if ~isempty(nbr_cells{j})
                    [min_c, min_index] = min(abs(cell_log_intensity(nbr_cells{j})-region_log_intensity(j)));
                    min_c_pairs(j, 1) = min_c;
                    min_c_pairs(j, 2) = nbr_cells{j}(min_index);
                else
                    min_c_pairs(j, 1) = realmax;
                    min_c_pairs(j, 2) = 0;
                end
            end
        end
    end
    
end
        
end