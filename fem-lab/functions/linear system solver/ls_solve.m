function [x] = ls_solve(A,b)
    
    [~,p] = chol(A)
    p=0
    if p==0
        %x = ls_cg(A,b,b);
    else
        %[x, error, iter, flag] = ls_gmres( A, 0*b, b, ichol(A), 100, 100, 10^-6 )
    end
    x=pcg(A,b,10^(-10),10000);
    ls_error=(A*x-b)'*(A*x-b)
end