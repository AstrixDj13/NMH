function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds = struct();

bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bnorth = {'PER'};
bconds.bsouth = {'PER'};

end

