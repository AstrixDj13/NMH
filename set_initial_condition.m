function [ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb )
%SET_INITIAL_CONDITION Set initial fields

% TODO: Bottom elevation
zb = -0.001 * grid.x' * ones(size(grid.y));

% TODO: Strickler coefficient
kst = 30 *ones( grid.M+2, grid.N+2 );

% TODO: Initial flow depth
h = 1 * ones( grid.M+2, grid.N+2 );

% TODO: Initial discharges
hu = kst*(sqrt(0.001)).*(h.^(5/3)) .* ones(grid.M+2, grid.N+2);
hv = 0 * ones(grid.M+2, grid.N+2);

end

