function x = ls_fgmres(A,b,m,tol,max_it,P,x0)    % We want to implement a flexible GMRES algorithm proposed by Youssef Saad.     % Input:    % m: restart length    % A: matrix of the linear system Ax=b    % b: right hand side of our linear system    % tol: tolerance for residual    % P: preconditioner    % x0: initial guess        % Output:    % x approximated solution of Ax=b        %Reference and Acknowledgements:    % 1) Saad: "A flexible inner-outer preconditioned GMRES algorithm", SIAM 1993 (Page 3: Algorithm 2.2 )    % 2) Nick Alger: http://scicomp.stackexchange.com/        % ------------------------- FMGRES ----------------------    tic    dim = size(b,1);    it = 1;    while (it <= max_it)        % Step 1) Start        H = zeros(m+1,m);                % Step 2) Arnoldi process        r0 = b - A*x0;        beta = norm(r0);        V = zeros(dim, m+1);        V(:,1) = (1/beta)*r0;        Z = zeros(dim,m);        for j=1:m            Z(:,j) = P*V(:,j);            w = A*Z(:,j);            for i=1:j                H(i,j) = w'*V(:,i);                w = w - H(i,j)*V(:,i);            endfor            H(j+1, j) = norm(w);            V(:, j+1) = (1/H(j+1, j))*w;                        % Step 3) Form the approximation solution            e1 = zeros(j+1,1);             e1(1) = 1;            y = H(1 : j+1, 1:j)\(beta*e1);            x = x0 + Z(:, 1:j)*y;            r_norm = norm(b-A*x);            if (r_norm < tol*norm(b))                disp(['FGMRES converged with relative tolerance ', ...                    num2str(r_norm/norm(b)), ...                    ' at iteration number ', ...                    num2str(it)])                return            endif            it = it + 1;            if (it > max_it)                break            endif        endfor            x0 = x;        % Step 4) Restart        endwhile    tocendfunction