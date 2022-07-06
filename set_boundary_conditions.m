function [ bconds ] = set_boundary_conditions()
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds = struct();

bconds.bwest = {'HFIX','HUFIX'};
bconds.beast = {'WALL'};
bconds.bnorth = {'PER'};
bconds.bsouth = {'PER'};

% TODO: Insert values for boundary condition
bconds.hwest = 1
bconds.huwest = 30*0.001^(1/2)
end

