function [ grid ] = generate_grid( grid )

%GENERATE_GRID Creates a struct containing all grid information.

% Compute grid spacing
grid.dx = (grid.xmax - grid.xmin) / grid.nx;
grid.dy = (grid.ymax - grid.ymin) / grid.ny;

% Create coordinate arrays in x-direction
grid.x(1) = grid.xmin - grid.dx/2;
grid.x(grid.nx+2) = grid.xmax + grid.dx/2;

for i = 1:grid.nx
    grid.x(i+1) = grid.x(i) + grid.dx;
end

% Create coordinate arrays in x-direction
grid.y(1) = grid.ymin - grid.dy/2;
grid.y(grid.ny+2) = grid.ymax + grid.dy/2;

for i = 1:grid.ny
    grid.y(i+1) = grid.y(i) + grid.dy;
end

