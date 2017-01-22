function [x] = ls_solve(A,b)

    % Calculate initial guess with pcg
    tic
    x= pcg(A,b,10^(-10),100);
    
    % Use initial guess for FGMRES
    x1= ls_fgmres(A,b,10,10^(-20),500,eye(rows(A)),x);
    gmres=toc
    tic
    x2=ls_cg(A,b,b);
    cg=toc
    % Error by linear solver
    ls_error=(A*x-b)'*(A*x-b)
    
end