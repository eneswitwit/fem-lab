function SM  = sm_assemble_global(mesh_size,SF,SM_local)    % This function is where we want to build the global stiffness matrix. For our approach that means that we have to    % compute the local stiffness matrix first and then build the global stiffness matrix via an indices renumbering like we     % saw in the lecture. For more information see the documentation.    % Initialize properties needed for renumbering the shape function and for the loops.    pol_deg=sqrt(rows(SF))-1;    cells_per_row= (1/mesh_size);    nodes_per_edge=pol_deg+1;    number_of_cells = cells_per_row^2;    number_of_nodes = ((pol_deg * cells_per_row)+1)^2;    SM = zeros(number_of_nodes,number_of_nodes);    SM2 = zeros(number_of_nodes,number_of_nodes);    Cell_nodes = mesh_cell(mesh_size,pol_deg);        % Compute local stiffness matrix    size_SM_local = length(SM_local);        % Build SM out of the entries of the local stiffness matrix via appropriate renumbering.    for k=1:number_of_cells        SM(Cell_nodes(k,:),Cell_nodes(k,:))+=SM_local;    endfor    endfunction