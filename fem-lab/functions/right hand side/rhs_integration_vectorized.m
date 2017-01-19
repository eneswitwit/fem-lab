function rhs = rhs_integration_vectorized(mesh_size,SF,f)
    tic
    %Let cell be the matrix, which stores the vertices for each cell
    %Let vertex be the matrix, which stores the coordinates for each vertex
    %Let SF be the matrix, containing the coefficients of the shape funtions
    
    %This function will give us the right hand side of the linear system
    
    % Useful computations for later use
    [vertex,cell]=mesh_generate(mesh_size);
    polynomial_deg = sqrt(length(SF))-1;
    nodes_per_cell=(polynomial_deg+1)^2;
    cells_per_row=(1/mesh_size);
    number_of_nodes=((polynomial_deg*cells_per_row)+1)^2;
    cell_nodes=mesh_cell(mesh_size,polynomial_deg);
    
    % Initialize Gauss Quadratur
    [sample_points,weights] = int_gauss_weights(polynomial_deg*10,0,1);
    
    rhs=zeros(number_of_nodes,1);
    
    for k=1:cells_per_row^2
    % This transformation is further explained in our documentation
        g =  @(x,y) f(mesh_size*x+vertex(cell(k,1),1),mesh_size*y+vertex(cell(k,1),2)).*hf_eval_poly(x,y,SF);
        integral=int_gauss(sample_points,weights,sample_points,weights,g);
        rhs(cell_nodes(k,:)')=rhs(cell_nodes(k,:)')+integral;
    endfor
    rhs=rhs*mesh_size^2;
    toc
endfunction