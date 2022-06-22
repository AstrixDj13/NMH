function [parm, flow] = set_params(parm, flow, infilename)

    % reading parameters from input file and parsing to strucs
    
    load(infilename, 'flowtype', 'dt', 'ntst', 'm', 'n', 'xl', 'yl', ...
        'nu', 'transporting_u', 'transporting_v');
            
    parm.flowtype = flowtype;
    parm.dt = dt;       
    parm.ntst = ntst;
    parm.m = m;         
    parm.n = n;
    parm.xl = xl;       
    parm.yl = yl;
    parm.nu = nu;
    flow.transporting_u = transporting_u;
    flow.transporting_v = transporting_v;
     
    
    %**********************************************************************
    % to manually override automatic input, uncomment below
    % and reset params:
    %
    % parm.dt = ...
    % parm.ntst = ...
    %
    %**********************************************************************
        
    
end




