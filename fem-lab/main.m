function [SM,rhs] = main()
    % Add all subfolders to working directory
    addpath(genpath([pwd '/functions']))
    % Initialize parameters mesh-size and polynomial degree
    mesh_size=1/2
    pol_deg=2
    % Initialize right hand side of strong formulation
    f = @(x,y) cos(x*pi).*cos(y*pi);
    %f=@(x,y) 3;
    % Initialize exact solution
    u_exact=@(x,y) cos(x*pi).*cos(y*pi)*(1/(1+2*pi*pi));
    
    [Vertex,Cell]=mesh_generate(mesh_size);
    SF=sf_generate(pol_deg);
    
    SM=sm_assemble_global(mesh_size,pol_deg);
    disp("------------------Assembled global matrix------------------");
    rhs=rhs_integration(mesh_size,SF,f);
    disp("------------------Assembled right hand side------------------");
    u_coeff=ls_solve(SM,rhs);
    disp("------------------Solved linear system------------------");
    m_plot_solution(u_coeff,pol_deg,mesh_size)
    %u=@(x,y) hf_eval_solution(x,y,u_coeff,Cell,Vertex,pol_deg,SF);
    %error_L2(u,u_exact,100)
    
    
    
endfunction