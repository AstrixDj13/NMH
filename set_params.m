function [grid, run, constants] = set_params(infilename)

    % Gravity [m/s^2]
    constants.g = 9.81;
    
    % number of ghost cell
    grid.NGHOST = 1;

    % reading parameters from input file and parsing to strucs

    load(infilename);
    
    grid.nx = nx;         
    grid.ny = ny;
    grid.xmax = xmax;       
    grid.xmin = xmin;
    grid.ymax = ymax;       
    grid.ymin = ymin;         
    
    run.dt = dt;       
    run.ntst = ntst;
     
end




