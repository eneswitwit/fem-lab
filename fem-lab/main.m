function main()
    % Add all subfolders to working directory
    genpath([pwd '\functions'],[pwd '\functions\bin'])
    addpath(genpath([pwd '\functions'],[pwd '\functions\bin']))
    % Initialize parameters mesh-size and polynomial degree
    mesh_size=0.5
    pol_deg=1
    % Initialize right hand side of strong formulation
    f = @(x,y) cos(x*pi).*cos(y*pi)
    
    [Vertex,Cell]=mesh_generate(mesh_size);
    SF=sf_generate(pol_deg);
    SM=sm_assemble_global(mesh_size,pol_deg);
    rhs=rhs_integration(mesh_size,pol_deg,Cell,Vertex,SF,f)
    
    
endfunction