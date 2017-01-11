function [x0] = ls_cg(A,b,x0)
    
    count=1;
    
    r0=b-A*x0;
    p0=r0;
    
    while (r0'*r0>10^(-6))
        alpha0=(r0'*r0)/(p0'*A*p0);
        %step x0 -> x1
        x0=x0+alpha0*p0;
        r1=r0-alpha0*A*p0;
        beta0=(r1'*r1)/(r0'*r0);
        %step r0 -> r1
        p0=r1 + beta0*p0;
        r0=r1;
        count++;
    endwhile 
    
endfunction
