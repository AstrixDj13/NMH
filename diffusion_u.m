function [diff_u] = diffusion_u(parm, flow, i, j)

    % returns the diffusive term for u at {x(i),y(j)}, i.e. the Laplacian
    % of u multiplied by the viscosity, using a central difference scheme
    % (CDS)
    
    % for convenience, we copy the following variables
    
    ip = parm.ip(i);
    im = parm.im(i);
    jp = parm.jp(j);
    jm = parm.jm(j);
    dx = parm.dx;
    dy = parm.dy;
    nu = parm.nu;
    
    % Now evaluate the diffusive term by applying CDS
    diff_u = (nu * (flow.u(ip,j) - 2*flow.u(i,j) + flow.u(im,j))/(parm.dx*parm.dx)) + (nu * (flow.u(i,jp) - 2*flow.u(i,j) + flow.u(i,jm))/(parm.dy*parm.dy));       
             
end