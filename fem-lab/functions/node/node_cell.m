% uWird nicht mehr gebrauchtfunction cell_matrix_node = node_cell(mesh_size,polynomial_deg)    dimension_domain = 2;    k = 0;    % nodes_counter gives us the number of nodes we except for a specific polynomial degree    nodes_per_cell =  (factorial(dimension_domain)/(factorial(k)*factorial(dimension_domain-k))) * polynomial_deg^k * (polynomial_deg+1)^(dimension_domain-k)    cell_counter = (1/mesh_size)^2    cell_matrix_node = zeros(cell_counter, nodes_per_cell);    row_counter = polynomial_deg + 1     nodes_per_row = nodes_per_cell / row_counter;    node_counter = 1;    node_per_cell_iteration=0;    k=1;    % mesh properties    mesh_cells_per_row=(1/mesh_size);    mesh_nodes_per_row=(polynomial_deg*mesh_cells_per_row)+1;    mesh_nodes_row_range = mesh_nodes_per_row;          for i=1:cell_counter        cell_matrix_node(i,1)=node_counter;        % check if at border        if(node_counter + nodes_per_row - 1 < mesh_nodes_row_range)            node_counter = node_counter + nodes_per_row - 1        elseif( i < cell_counter )            node_counter = node_counter + nodes_per_row - 1 - nodes_per_row*mesh_cells_per_row+mesh_cells_per_row+mesh_nodes_per_row*(row_counter-1);            mesh_nodes_row_range = mesh_nodes_row_range + mesh_nodes_per_row*(row_counter-1);        endif    endfor        for i=1:cell_counter        for k=2:nodes_per_row            cell_matrix_node(i,k)=cell_matrix_node(i,k-1)+1;        endfor    endfor        for c=1:cell_counter        for k=1:nodes_per_row            i=k+nodes_per_row;            while i <= nodes_per_cell                cell_matrix_node(c,i)=cell_matrix_node(c,i-nodes_per_row)+mesh_nodes_per_row;                i=i+nodes_per_row;            endwhile          endfor      endforendfunction 