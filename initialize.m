function [parm, flow] = initialize(parm, flow)

    % initialization of indices: ip, im, jp, jm, 
    % meaning respectively: i+1, i-1, j+1 and j-1
    
    % Note: the spatial discretization scheme and/or boundary condition 
    % can be controlled by these index fields on the boundaries
    
    
    for i = 1 : parm.m

        parm.ip(i) = i+1;
        parm.im(i) = i-1;

    end
    
    for j = 1 : parm.n

        parm.jp(j) = j+1;
        parm.jm(j) = j-1;

    end
    
    % set periodic boundary conditions
   
    parm.ip(parm.m) = 2;  
    parm.jp(parm.n) = 2;
    parm.im(1) = parm.m - 1;
    parm.jm(1) = parm.n - 1;
    
    % compute and set dx, dy
    
    parm.dx = parm.xl/(parm.m-1);
    parm.dy = parm.yl/(parm.n-1);
          
    % set initial condition for the velocity field
    [parm, flow] = initcond_2d_condiff_scal(parm, flow);
               
end

function [parm, flow] = initcond_2d_condiff_scal(parm, flow)

% define the initial velocity field u for the case of the 2D
% convection-diffusion equation with constant transporting velocities

    for i = 1 : parm.m
        for j = 1 : parm.n
            
           flow.u(i,j) = mod((i > floor(parm.m/2)) + (j > floor(parm.n/2)), 2);  
           
        end
    end
end

