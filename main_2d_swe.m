%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 7
%
% This code solves the 2D Shallow Water equations
%
% author: H. Zeng & L. Unglehrt
% June, 2020
%**************************************************************************
%% Set up directory for results
clear;
close all

if ~exist('results', 'dir')
       mkdir('results')
    else
        delete('results\*')
 end

%% Initialize simulation
% read infile 
infilename = 'infile_2D_swe_test.mat';
load(infilename)
fprintf('infilename is: %s\n', infilename)

% build structures 
[grid, run, constants, flow, bconds] = build_structs;
fprintf('struct built\n')

% fill some fields with data from input file
[grid, run, constants] = set_params(infilename);
fprintf('parameters set\n')

% Generate an equidistant grid 
[grid] = generate_grid(grid);    

% Set initial conditions 
run.t = 0;
[ flow ] = set_initial_condition( grid, flow );

figure(1)
surf(flow.h + flow.zb)
colorbar
xlabel("x")
ylabel("y")
zlabel("Water Level h")
print('-djpeg', 'results\initial_surface','-r250')

figure(2)
contourf(flow.h + flow.zb)
a = colorbar;
a.Label.String = "Water level h";
xlabel("x")
ylabel("y")
zlabel("Water Level h")
print('-djpeg', 'results\initial_contour','-r250')

% Create boundary conditions
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};

waterLevel = zeros(length(flow.h), length(flow.h));


 

%% Time integration
count=10;

for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );

% Plot results after 10th time step
 if mod(itstep,10) == 0
    height_max(itstep/10) = max(max(flow.h + flow.zb))
    height_min(itstep/10) = min(min(flow.h + flow.zb))
    
    figure(3)
    surf(flow.h + flow.zb)
    colorbar
    caxis([0.9 1])
    xlabel("x")
    ylabel("y")
    zlabel("Water Level (h)")
    zlim([0.0 2.0])
    print('-djpeg', 'results\timestep' + string(count),'-r250')
    count= count+10;

 end

end

figure(5)
plot(height_max)
hold on
plot(height_min)
xlabel("time step")
ylabel("min/max water level h of whole domain")
hold off
legend({'max h','min h'},'Location','northeast')
print('-djpeg', 'results\min_max_waterlevel','-r250')