function [ bconds ] = set_boundary_conditions(constants)
%SET_BOUNDARY_CONDITIONS Set boundary conditions

bconds = struct();

bconds.bwest = {'HFIX','HUFIX'};
bconds.beast = {'WALL'};
bconds.bnorth = {'PER'};
bconds.bsouth = {'PER'};

% TODO: Insert values for boundary condition
bconds.hwest = constants.h
bconds.huwest = constants.kst*constants.I^(1/2)*constants.h^(5/3)
end

