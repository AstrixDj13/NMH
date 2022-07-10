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

% Set up directory for results
clear;
close all

if ~exist('results', 'dir')
       mkdir('results')
       mkdir("results\Flowdepth")
       mkdir("results\Waterlevel")
       mkdir("results\Discharge")
       mkdir("results\Waveheight")
       mkdir("results\Wavespeed")
    else
        delete('results\Discharge\*')
        delete('results\Flowdepth\*')
        delete('results\Waterlevel\*')
        delete('results\Waveheight\*')
        delete('results\Wavespeed\*')
 end

if is_octave()
    % ---- Load packages ----------------------------------------------------
    pkg load statistics
end

% ---- Constants ----------------------------------------------------------
constants = struct();
constants.g = 9.81;
constants.enable_entropy_fix = true;
constants.NGHOST = 1;
constants.I = 0.001;
constants.kst = 30;
constants.h = 1;
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
[ h, hu, hv, kst, zb ] = set_initial_condition( grid, h, hu, hv, kst, zb, constants );

% ---- Create boundary conditions -----------------------------------------
bconds = set_boundary_conditions(constants);

% ---- Computation --------------------------------------------------------
% TODO: Number of time steps
mtstep = 200;

% TODO: Time step size
dt = 0.01;

% Frequency of diagnostic output
itdiag = 1;

% Initial location and time of wave
count =10;
x(1) = 10;
t_a(1) = 0;
a(1)=0;
time = 0;

%% Time integration
for itstep = 1:mtstep
    % Reset time stepping function (because of persistent variables)
    [ t, h, hu, hv ] = time_step_rk(itstep==1, constants, grid, dt, ...
        t, h, hu, hv, kst, zb, bconds );
    
    % Diagonostic output
    if mod(itstep, itdiag) == 0
        % CFL number
        [CFL, FR] = compute_CFL_number(constants, grid, dt, h, hu, hv);
        fprintf('%d : CFL number:         %e     Froude number:         %e\n', itstep, ...
            CFL, FR);
    end
    
    
    % Water surface
    zw = zb + h;
    
    % Wave height
    height(itstep)= max(zw(:,2))-1;
    % Midpoint
    midpoint(itstep) = max(zw(:,2))-height(itstep)/2 ;
    % Index of midpoint
    k = find(zw(:,2)>=midpoint(itstep));
    % Position on x-Axis of midpoint
    x(itstep+1)= (k(1)-1)*grid.dx;
    % Time vector
    t_a(itstep+1)= itstep*dt;
    % If change in position save time and compute velocity
    if x(itstep+1)~=x(itstep)
        dtime = t_a(itstep+1) - time;
        time = t_a(itstep+1);
        a(itstep+1)=(x(itstep+1)-x(itstep))/dtime;
    else
        a(itstep+1)=a(itstep);
     
    end
    

    
    f1=figure(1);
    plot(grid.x, h(:,2))
    xlabel('x')
    ylabel('h')
    title(sprintf('Flow depth at t=%g',t))
    xlim([0,10])
    % ylim([0,4])
    drawnow
    
    f2=figure(2);
    plot(grid.x, zw(:,2))
    xlabel('x')
    ylabel('zw')
    title(sprintf('Water level at t=%g, CFL=%.2f',t,CFL))
    xlim([0,10])
    % ylim([0,4])
    drawnow
    
    f3=figure(3);
    plot(grid.x, hu(:,2), 'r')
    xlabel('x')
    ylabel('hu')
    title(sprintf('Specific discharge at t=%g',t))
    xlim([0,10])
    %ylim([0,5])
    drawnow

    if mod(itstep,10) == 0
        print(f2,'-djpeg', 'results\Waterlevel\Waterlevel_timestep' + string(count),'-r250')
        print(f1,'-djpeg', 'results\Flowdepth\Flowdepth_timestep' + string(count),'-r250')
        print(f3,'-djpeg', 'results\Discharge\Discharge_timestep' + string(count),'-r250')
        count = count+10;
    end
end

figure(4)
plot(height)
xlabel('timestep')
ylabel('Wave height')
print('-djpeg', 'results\Waveheight\Waveheight','-r250')

figure(5)
plot(a)
xlabel('timestep')
ylabel('Wave speed [m/s]')
xlim([0 200])
print('-djpeg', 'results\Wavespeed\Wavespeed','-r250')
