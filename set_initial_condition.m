function [ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb )
%SET_INITIAL_CONDITION Set initial fields

% TODO: Bottom elevation
 zb(:,:) = 0.01 * grid.x' * ones(size(grid.y));

% TODO: Strickler coefficient
 kst(:,:) = 30;

% TODO: Initial flow depth
 h(:,:) = ((1*(grid.x'<=2)) + (0.1*(grid.x'>2)))* ones(size(grid.y));

% TODO: Initial discharges
 hu(:,:) = 0;
 hv(:,:) = 0;
end

