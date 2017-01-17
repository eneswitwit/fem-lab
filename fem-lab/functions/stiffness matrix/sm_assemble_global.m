function SM  = sm_assemble_global(mesh_size,polynomial_deg)    % Initialize properties needed for renumbering the shape function and for the loops    cells_per_row= (1/mesh_size);    nodes_per_edge=polynomial_deg+1;    number_of_cells = cells_per_row^2;    number_of_nodes = ((polynomial_deg * cells_per_row)+1)^2;    SF = sf_generate(polynomial_deg);    SM = zeros(number_of_nodes,number_of_nodes);    cell_matrix = mesh_cell(mesh_size,polynomial_deg);        % Calculate local stiffness matrix    SM_local = sm_assemble_local(SF,mesh_size);    size_SM_local = length(SM_local);        % Build SM out of the entries of the local stiffness matrix via appropriate renumbering    for k=1:number_of_cells        for i=1:nodes_per_edge^2            for j=1:nodes_per_edge^2                SM(cell_matrix(k,i),cell_matrix(k,j))=SM(cell_matrix(k,i),cell_matrix(k,j))+SM_local(i,j);            endfor        endfor    endforendfunction