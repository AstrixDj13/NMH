function [flow] = euler_2d_condiff_scal(parm, flow)

    % Perform one time integration step by the Explicit Euler Scheme, viz.
    % phi^n+1 = phi^n + dt * rhs(phi^n)

    % Tip: MATLAB syntax does not require loops in order to perform elementwise 
    % matrix operations, so the implemantation is simple ...

    flow.u = flow.u + parm.dt * flow.rhsu;
    
end

