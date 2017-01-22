function main()
    tic
    % Add all subfolders to working directory
    addpath(genpath([pwd '\functions']))
    % Initialize parameters mesh-size and polynomial degree
    mesh_size=1/64;
    pol_deg=1;
    % Initialize right hand side of strong formulation
    f = @(x,y) cos(x*pi).*cos(y*pi);
    % Initialize exact solution
    u_exact=@(x,y) cos(x*pi).*cos(y*pi)*(1/(1+2*pi*pi));
    
    % Initialize our mesh and our coefficient matrix for the shape functions
    [Vertex,Cell]=mesh_generate(mesh_size);
    SF=sf_generate(pol_deg);
    
    % sm_assemble_global will give the global stiffness matrix
    % sm_assemble_global will call the function sm_assemble_local, which initializes the local stiffness matrix
    SM=sm_assemble_global(mesh_size,pol_deg);
    disp("------------------Assembled global matrix------------------");
    % Initialize right hand side of our linear system
    rhs=rhs_integration(mesh_size,SF,f);
    disp("------------------Assembled right hand side------------------");
    % solve the linear system using ls_solve
    u_coeff=ls_solve(SM,rhs);
    disp("------------------Solved linear system------------------");
    % hf_eval_solution is used to evaluate the approximation, given the coefficients u_coeff;
    u=@(x,y) hf_eval_solution(x,y,u_coeff,Cell,Vertex,pol_deg,SF);
    % m_plot will plot the approximation, the exact solution and the right hand side of the strong form to illustrate, that it is in fact a multiple of the solution
    m_plot_solution(u_coeff,pol_deg,mesh_size)
    % error_L2 computes the error, using the exact solution and the L2-norm
    error_L2(u,u_exact,10)
    % overall runtime using 'tic' 'toc'
    overall_runtime=toc
    
    
endfunction