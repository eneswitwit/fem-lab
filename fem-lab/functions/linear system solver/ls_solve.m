function [x] = ls_solve(A,b)

    % Calculate initial guess with pcg
    x= pcg(A,b,10^(-10),100);
    
    % Use initial guess for FGMRES
    x= ls_fgmres(A,b,10,10^(-20),500,eye(rows(A)),x);
    
    % Error by linear solver
    ls_error=(A*x-b)'*(A*x-b)
    
end