function [ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb, constants )
%SET_INITIAL_CONDITION Set initial fields

% TODO: Bottom elevation
zb = -constants.I * grid.x' * ones(size(grid.y));

% TODO: Strickler coefficient
kst = constants.kst *ones( grid.M+2, grid.N+2 );

% TODO: Initial flow depth
h = constants.h * ones( grid.M+2, grid.N+2 );

% TODO: Initial discharges
hu = kst*(sqrt(constants.I)).*(h.^(5/3)) .* ones(grid.M+2, grid.N+2);
hv = 0 * ones(grid.M+2, grid.N+2);

end
