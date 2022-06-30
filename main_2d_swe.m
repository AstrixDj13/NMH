%**************************************************************************
% NMH Lab Summer Semester 2020, Assignment 7
%
% This code solves the 2D Shallow Water equations
%
% author: H. Zeng & L. Unglehrt
% June, 2020
%**************************************************************************
clear;
close all

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
zlabel("Water Level (h)")
% Create boundary conditions
bconds.bwest = {'WALL'};
bconds.beast = {'WALL'};
bconds.bsouth = {'WALL'};
bconds.bnorth = {'WALL'};


%% Time integration
for itstep = 1:run.ntst
    [ run, flow ] = time_step_rk( itstep==1, constants, grid, run, ...
        flow, bconds );
 %for itstep = 1:10:run.ntst
 if mod(itstep,10) == 0

    height_max(itstep/10) = max(max(flow.h + flow.zb))
    height_min(itstep/10) = min(min(flow.h + flow.zb))
 end
%% Plot results
% TODO TODO TODO TODO TODO TODO TODO
% figure(2)
% surf(flow.h + flow.zb)
% colorbar
% xlabel("x")
% ylabel("y")
% zlabel("Water Level (h)")

end
figure(2)
plot(height_max)
hold on
plot(height_min)
hold on
plot(height_max - height_min)
hold off