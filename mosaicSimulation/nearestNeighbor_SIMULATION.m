function [nnDistances] = nearestNeighbor_SIMULATION(centers)
    nnDistances = zeros(size(centers,1), 1);
    for i = 1:size(centers, 1)
        center_i = centers(i, :);
        minDistance = 0;
        for j = 1:size(centers, 1)
            if j == i
                continue
            end

            comp_j = centers(j, :);
            distance_j = pdist([center_i; comp_j], 'euclidean');
            if minDistance == 0
                minDistance = distance_j;
            elseif minDistance > distance_j
                minDistance = distance_j;
            end
        end
        nnDistances(i) = minDistance;
    end

end