function [parm, flow] = timeint(parm, flow)
   
    % loop over timesteps 
    for itst = 1 : parm.ntst

        disp(['time step ', num2str(itst)]);

        % compute right hand side
        [flow.rhsu] = rhs_2d_condiff_scal(parm, flow);

        % perform one explicit Euler time step
        [flow] = euler_2d_condiff_scal(parm, flow);

    end
            
end

