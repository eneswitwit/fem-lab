function [x] = ls_solve(A,b)
    tic
    x = ls_cg(A,b,0*b);
    toc
end