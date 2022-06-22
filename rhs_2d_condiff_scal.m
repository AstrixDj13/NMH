function [rhsu] = rhs_2d_condiff_scal(parm, flow)

    % local initialization
    
    rhsu = zeros(parm.m, parm.n);
    
    % computation of right hand side in a loop for every grid point
    % Note: diffusion_u(...) is in a separate file 

    for i = 1 : parm.m
        for j = 1 : parm.n
            
            % compute the convective (advection) term
            ad = advection_u_lin2d(parm, flow, i, j);
            
            % compute the diffusion term
            diff = diffusion_u(parm, flow, i, j);
            
            rhsu(i,j) = ad + diff;
                
        end
    end

end

% used functions

function [adv_u] = advection_u_lin2d(parm, flow, i, j)

    % returns the advective term of u at {x(i),y(j)}, using a 
    % central difference scheme (CDS)
    
    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    
    adv_u = -  flow.transporting_u*(flow.u(ip,j) - flow.u(im,j))/(2*parm.dx)  - flow.transporting_v*(flow.u(i,jp) - flow.u(i,jm))/(2*parm.dy);
        
end