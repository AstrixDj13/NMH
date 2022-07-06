function [CFL_max] = compute_CFL_number(constants, grid, dt, h, hu, hv)
%COMPUTE_CFL_NUMBER Computes the CFL number. This is probably not perfect.


CFL_max = 0;

for j = grid.NGHOST+1:length(grid.y)-grid.NGHOST
    % The i-loop is inside the j-loop because of the column major memory
    % layour in MATLAB
    for i = grid.NGHOST+1:length(grid.x)-grid.NGHOST
        % compute velocity at cell center
        u = hu(i,j) / h(i,j);
        v = hv(i,j) / h(i,j);
        c = (constants.g*h(i,j))^(0.5);
        % TODO: Compute CFL
        u_res = max(abs(u+c),abs(u-c));
        CFL = (u_res*dt/grid.dx) + (v*dt/grid.dy);

        % use worst case
        CFL_max = max( CFL, CFL_max );

    end
end

end

