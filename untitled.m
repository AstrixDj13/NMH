
    for i = 1 : parm.m
        for j = 1 : parm.n
            
           flow.u(i,j) = mod((i > floor(parm.m/2)) + (j > floor(parm.n/2)), 2);  
           
        end
    end