function linear_solver()        for i=1:100        % Calculate Stiffness Matrix        mesh_size= 1/i ;        pol_deg= 1;        % Initialize our mesh and our coefficient matrix for the shape functions        [Vertex,Cell]=mesh_generate(mesh_size);        SF=sf_generate(pol_deg);        % sm_assemble_local computes the local stiffness matrix. In our case,the local stiffness matrix looks the same for every cell.        SM_local=sm_assemble_local(mesh_size,SF);        % sm_assemble_global will give the global stiffness matrix        SM=sm_assemble_global(mesh_size,SF,SM_local);            % FGMRES        tic        % Preconditioning Matrix: Incomplete Cholesky factorization with no fill-in        M = ichol(SM);        x1= ls_fgmres(A,b,100,10^(-10),50,M,b);        gmres_runtime=toc;            % MINRES        tic        x2= ls_minres(A,b,b,100,10^(-10))        minres_runtime=toc;            % CG Linear solver        tic        x3=ls_cg(A,b,b);        cg_runtime=toc;            % Error by linear solver        ls_error_fgm=(A*x1-b)'*(A*x1-b);        ls_error_cg=(A*x3-b)'*(A*x3-b);        ls_error_minres=(A*x2-b)'*(A*x2-b);            endfor        end