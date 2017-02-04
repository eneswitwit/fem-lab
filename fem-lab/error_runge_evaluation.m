function error_runge_evaluation()    add_fem_lab();    % --- Acknowledgments ---    % 'A posteriori analysis for PDE's' by S. Repin, Page 60    % Document can be found here: https://math.aalto.fi/opetus/lme/ApostPDERepin.pdf    % -----------------------        % The Goal of that function is to show that the old fashioned error estimator from C. Runge is actually a good error estimator for our    % case. Therefore we want to have a sequence of mesh sizes like    % (mesh_size, mesh_size/2, mesh_size/4, mesh_size/8, ...)     % and show that the solution is converging with decreasing mesh size.    % So we want to show computationally that the sequence of solutions with decreasing mesh sizes is a cauchy sequence.    % Also we want to show that the convergence order is similar to the convergence order of the L2 error with the analytical solution    % converging to 0 with decreasing mesh size.    % Initialize necessary variables, constants and functions    MIN_MESH_SIZE = 48;    counter = 1;    pol_deg= 1;    SF=sf_generate(pol_deg);        % Initialize vectors containg the errors of L2 and the error estimator with Runge     l2 = zeros(MIN_MESH_SIZE,1);    runge = zeros(MIN_MESH_SIZE,1);                % Initialize right hand side of strong formulation    f = @(x,y) cos(x*pi).*cos(y*pi);    % Initialize exact solution    u_exact= @(x,y) cos(x*pi).*cos(y*pi)*(1/(1+2*pi*pi));    % Initial guess of the solution is constant 1.    u_1 = @(x,y) 1;    mesh_size = 1;        while(mesh_size > 1/MIN_MESH_SIZE )                % Mesh size for u_1 is equal to mesh_size and mesh_size for u_2 should be mesh_size/2         mesh_size = mesh_size/2;                % Initialize vertex,cell and shape function matrices for evaluating the solution and solving FEM        [Vertex_2,Cell_2]=mesh_generate(mesh_size/2);        SF = sf_generate(pol_deg);            % Compute solution for mesh_size/2        SM_local=sm_assemble_local_vectorized(mesh_size/2,SF);        SM =sm_assemble_global(mesh_size/2,pol_deg,SM_local);        rhs=rhs_integration_vectorized(Vertex_2,Cell_2,SF,f);        u_coeff_2 = ls_solve(SM,rhs);        u_2=@(x,y) hf_eval_solution(x,y,u_coeff_2,Cell_2,Vertex_2,pol_deg,SF);                % Compute L2 error of solution with mesh size equal to mesh_size and the exact solution        l2(counter) = error_L2( u_1 , u_exact , 50)                % Compute the L2 error of the solution with mesh size equals to mesh_size and the solution with mesh size equals to mesh_size/2        runge(counter) = error_L2(u_1,u_2,50)                counter++;        u_1 = u_2;            endwhile        % Plot results    x = zeros(1,size(l2));    mesh_size = MAX_MESH;    for i=1:size(l2)        x(i) = mesh_size;        mesh_size = mesh_size/2;    endfor    plot(x,l2,'r',x,runge,'b');endfunction   