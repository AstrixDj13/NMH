function [ grid ] = generate_grid( xmin, xmax, ymin, ymax, nx, ny, NGHOST )
%GENERATE_GRID Creates a struct containing all grid information.

% Compute grid spacing
dx = ( xmax - xmin ) / nx;
dy = ( ymax - ymin ) / ny;

% Create coordinate arrays
x = ( xmin - ( NGHOST - 0.5 ) * dx ):dx:( xmax + ( NGHOST - 0.5 ) * dx );
y = ( ymin - ( NGHOST - 0.5 ) * dy ):dy:( ymax + ( NGHOST - 0.5 ) * dy );

% Create empty struct
grid = struct();

grid.dx = dx;
grid.dy = dy;
grid.x = x;
grid.y = y;
grid.NGHOST = int32(NGHOST);
grid.M = nx;
grid.N = ny;

end

