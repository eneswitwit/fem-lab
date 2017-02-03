function rhs = rhs_integration_new(Vertex,Cell,SF,f)

    tic

    %Let cell be the matrix, which stores the vertices for each cell
    %Let vertex be the matrix, which stores the coordinates for each vertex
    %Let SF be the matrix, containing the coefficients of the shape funtions
    
    %This function will give us the right hand side of the linear system
    
    % Initialize Gauss Quadratur
    [sample_points,weights] = int_gauss_weights(10,0,1);
    % Useful computations for later use
    mesh_size=Vertex(2,1)-Vertex(1,1);
    pol_deg = sqrt(length(SF))-1;
    cells_per_row=(1/mesh_size);
    number_of_cells=cells_per_row^2;
    number_of_nodes=((pol_deg*cells_per_row)+1)^2;
    
    Nodes=mesh_cell(mesh_size,pol_deg);
    
    rhs=zeros(number_of_nodes,1);
    
    for k=1:number_of_cells
        active_vertex=Vertex(Cell(k,1),:);
        g =  @(x,y) repmat(f(mesh_size*x+active_vertex(1),mesh_size*y+active_vertex(2)),length(SF),1).*hf_eval_poly(x,y,SF);
        rhs(Nodes(k,:))+=int_gauss_vectorized_matrices(sample_points,weights,sample_points,weights,g);
    endfor
    
    rhs*=(mesh_size^2);
    
    new=toc
    
endfunction