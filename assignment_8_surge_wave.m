%**************************************************************************
% NMH Lab Summer Semester 2022
%
% Assignment 8: Surge wave
%
% author: L. Unglehrt, Y. Sakai
% 04.07.2022
%**************************************************************************
clear
close all
clc

if is_octave()
    % ---- Load packages ----------------------------------------------------
    pkg load statistics
end

% ---- Constants ----------------------------------------------------------
constants = struct();
constants.g = 9.81;
constants.enable_entropy_fix = true;
constants.NGHOST = 1;

% ---- Create initial conditions ------------------------------------------
t = 0.0;

% Number of cells
M = 125;
N = 1;

% Grid (xmin, xmax, ymin, ymax, nx, ny, nghost)
grid = generate_grid(0, 10, -1, 1, M, N, constants.NGHOST);

% Create empty fields
h   = zeros( length(grid.x), length(grid.y) );
hu  = zeros( length(grid.x), length(grid.y) );
hv  = zeros( length(grid.x), length(grid.y) );
zb  = zeros( length(grid.x), length(grid.y) );
kst = zeros( length(grid.x), length(grid.y) );

% Set initial condition
[ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb );

% ---- Create boundary conditions -----------------------------------------
bconds = set_boundary_conditions();

% ---- Computation --------------------------------------------------------
% TODO: Number of time steps
mtstep = 200;

% TODO: Time step size
dt = 0.01;

% Frequency of diagnostic output
itdiag = 1;

%% Time integration
for itstep = 1:mtstep
    % Reset time stepping function (because of persistent variables)
    [ t, h, hu, hv ] = time_step_rk(itstep==1, constants, grid, dt, ...
        t, h, hu, hv, kst, zb, bconds );
    
    % Diagonostic output
    if mod(itstep, itdiag) == 0
        % CFL number
        fprintf('%d : CFL number:         %e\n', itstep, ...
            compute_CFL_number(constants, grid, dt, h, hu, hv));
    end
    
    
    % Water surface
    zw = zb + h;
  
    figure(1)
    plot(grid.x, h(:,2))
    xlabel('x')
    ylabel('h')
    title(sprintf('Flow depth at t=%g',t))
    xlim([0,10])
    % ylim([0,4])
    drawnow
    
    figure(2)
    plot(grid.x, zw(:,2))
    xlabel('x')
    ylabel('h')
    title(sprintf('Water level at t=%g',t))
    xlim([0,10])
    % ylim([0,4])
    drawnow
    
    figure(3)
    plot(grid.x, hu(:,2), 'r')
    xlabel('x')
    ylabel('hu')
    title(sprintf('Specific discharge at t=%g',t))
    xlim([0,10])
    %ylim([0,5])
    drawnow
    
end
max = max(max(zw))
