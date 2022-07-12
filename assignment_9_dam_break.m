%**************************************************************************
% NMH Lab Summer Semester 2022
%
% Assignment 9: Dam break
%
% author: L. Unglehrt, Y. Sakai
% 11.07.2022
%**************************************************************************

% This test case checks whether the transition from subcritical to
% supercritical flow is represented correctly.

if is_octave()
    % ---- Load packages ----------------------------------------------------
    pkg load statistics
end

% ---- Constants ----------------------------------------------------------
constants = struct();
constants.g = 9.81;

% TODO
constants.enable_entropy_fix = true;

% ---- Create initial conditions ------------------------------------------
t = 0.0;

% Number of cells
M = 125;
N = 1;

% Grid (xmin, xmax, ymin, ymax, nx, ny, nghost)
grid = generate_grid(0, 10, -1, 1, M, N, 1);

% Create empty fields
h   = zeros( length(grid.x), length(grid.y) );
hu  = zeros( length(grid.x), length(grid.y) );
hv  = zeros( length(grid.x), length(grid.y) );
zb  = zeros( length(grid.x), length(grid.y) );
kst = zeros( length(grid.x), length(grid.y) );

% TODO: Set initial condition
[ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb );

% ---- Create boundary conditions -----------------------------------------
% TODO
bconds = set_boundary_conditions();

% ---- Compute initial water volume ---------------------------------------
initialVolume = ( grid.dx * grid.dy ) * sum( sum( h(2:end-1,2:end-1), 1 ), 2 );

% ---- Computation --------------------------------------------------------
% TODO: Number of time steps
mtstep = 400;

% Time step size
dt = 0.005;

% Frequency of diagnostic output
itdiag = 10;

%% Time integration
% Reset time stepping function (because of persistent variables)
for itstep = 1:mtstep
    [ t, h, hu, hv ] = time_step_rk(itstep==1, constants, grid, dt, ...
        t, h, hu, hv, kst, zb, bconds );

    % Diagonostic output
    if mod(itstep, itdiag) == 0
        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep, ...
            compute_CFL_number(constants, grid, dt, h, hu, hv));

        % Total volume
        fprintf('%d : Total fluid volume: %e\n', itstep, ...
            sum( sum( h(2:end-1,2:end-1), 1), 2 ) * grid.dx * grid.dy );

        % Mean energy head
        w = zb(2:end-1,2:end-1) + h(2:end-1,2:end-1);
        e = w + (hu(2:end-1,2:end-1).^2 + hv(2:end-1,2:end-1).^2 )./ ...
            (2 *constants.g * h(2:end-1,2:end-1).^2);
        fprintf('%d : Mean energy head: %e\n', itstep, ...
            mean( reshape( e, [], 1 ) ) );
        
       % Froude number
        u = hu(2:end-1,2)./h(2:end-1,2);
        Fr = u./sqrt(constants.g * h(2:end-1,2));
        

        % Water surface
        figure(1)
        plot(grid.x(2:end-1), h(2:end-1,2))
        xlabel('x')
        ylabel('h')
        title(sprintf('Flow depth at t=%g',t))
        ylim([0,1.1])
        drawnow

        % Specific discharge
        figure(2)
        plot(grid.x(2:end-1), hu(2:end-1,2), 'r')
        xlabel('x')
        ylabel('hu')
        title(sprintf('Specific discharge at t=%g',t))
        ylim([0,1.1])
        drawnow

        % TODO: Plot energy head, Froude number
        figure(3)
        plot(grid.x(2:end-1), e)
        xlabel('x')
        ylabel('Energy Head')
        title(sprintf('Energy Head at t=%g',t))
        ylim([0,1])
        drawnow

        figure(4)
        plot(grid.x(2:end-1), Fr)
        xlabel('x')
        ylabel('Fr')
        title(sprintf('Froude number at t=%g',t))
        %ylim([0,1])
        drawnow
    end

end

% ---- Compute final water volume -----------------------------------------
finalVolume = ( grid.dx * grid.dy ) * sum ( h(2:end-1,2:end-1), 'all' );
