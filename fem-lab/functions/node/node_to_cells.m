function node_belongs_to_cells_vector = node_to_cells(i_th_node,vertex_matrix,cell_matrix)    node_belongs_to_cells_vector = zeros(1,4);    %check if i-th node is a vertex in first column of cell-matrix    for i=1:length(cell_matrix)        if(i_th_node == cell_matrix(i,1))            node_belongs_to_cells_vector(1,1) = i;            break;        endif    endfor    %check if i-th node is a vertex in 2nd column of cell-matrix    for i=1:length(cell_matrix)        if(i_th_node == cell_matrix(i,2))            node_belongs_to_cells_vector(1,2) = i;            break;        endif    endfor    %check if i-th node is a vertex in 3rd column of cell-matrix    for i=1:length(cell_matrix)        if(i_th_node == cell_matrix(i,3))            node_belongs_to_cells_vector(1,3) = i;            break;        endif    endfor    %check if i-th node is a vertex in 4th column of cell-matrix    for i=1:length(cell_matrix)        if(i_th_node == cell_matrix(i,4))            node_belongs_to_cells_vector(1,4) = i;            break;        endif    endforendfunction    