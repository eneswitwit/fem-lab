function [x,ls_error_cg,ls_error_fgm,cg_runtime,gmres_runtime] = ls_solve(A,b)

    % FGMRES
    tic
    % Calculate initial guess with pcg
    x= pcg(A,b,10^(-10),100);
    % Use initial guess for FGMRES
    x1= ls_fgmres(A,b,10,10^(-10),50,eye(rows(A)),b);
    gmres_runtime=toc;
    
    % CG Linear solver
    tic
    x2=ls_cg(A,b,b);
    cg_runtime=toc;
    
    % Error by linear solver
    ls_error_fgm=(A*x1-b)'*(A*x1-b)
    ls_error_cg=(A*x2-b)'*(A*x2-b)
    
    if (ls_error_cg<ls_error_fgm)
        x=x2;
    else
        x=x1;
    endif

    
end