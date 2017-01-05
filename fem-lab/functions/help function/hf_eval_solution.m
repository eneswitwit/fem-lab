function val=hf_eval_solution(x,y,u,cell_matrix,vertex_matrix,polynomial_deg,SF)

    mesh_size=vertex_matrix(2,1)-vertex_matrix(1,1)
    cells_per_row=(1/mesh_size);
    nodes_per_row=(polynomial_deg*cells_per_row)+1;
    nodes_per_edge=polynomial_deg+1;
    
    val=0;
    
    if (floor(x*(nodes_per_row-1))==x*(nodes_per_row-1) && floor(y*(nodes_per_row-1))==y*(nodes_per_row-1))
        node_i=x*(nodes_per_row-1)
        node_j=y*(nodes_per_row-1)
        node_no=node_j*nodes_per_row+node_i+1
        val=u(node_no)
    else
        if (floor(x*cells_per_row)==x*cells_per_row)        % x on edge
            cell_j=floor(y*cells_per_row)
            if(x==0)                                                    % x on lower boundary
                cell_i=0
                active_cell=(cell_j)*cells_per_row+1
            else                                                   % x not on lower boundary
                if(x==1)                                                    % x on upper boundary
                    cell_i=cells_per_row-1
                    active_cell=(cell_j+1)*cells_per_row
                else                                                   % x not on lower boundary -> on an interior edge
                    cell_i=[x*cells_per_row-1;x*cells_per_row]
                    active_cell=[cell_j*cells_per_row+cell_i(1)+1;cell_j*cells_per_row+cell_i(2)+1]
                    cell_j=[cell_j;cell_j]
                endif
            endif
        else                                                   % x not on edge
            if (floor(y*cells_per_row)==y*cells_per_row)            % y on edge
                cell_i=floor(x*cells_per_row)
                if (y==0)                               % y on lower boundary
                    cell_j=0
                    active_cell=cell_i+1
                else                                                    % y not on lower boundary
                    if (y==1)                     % y on upper boundary
                        cell_j=cells_per_row-1
                        active_cell=cells_per_row*(cells_per_row-1)+cell_i+1
                    else                                                      % y not on lower boundary -> y on an interior edge
                        cell_j=[y*cells_per_row-1;y*cells_per_row]
                        active_cell=[cell_j(1)*cells_per_row+cell_i+1;cell_j(2)*cells_per_row+cell_i+1]
                        cell_i=[cell_i;cell_i]
                    endif
                endif
            else                                                    % y not on any edge -> (x,y) in the interior
                cell_i=floor(x*cells_per_row)
                cell_j=floor(y*cells_per_row)
                active_cell=cell_j*cells_per_row+cell_i+1
            endif
        endif
        active_cell
        for i=1:length(active_cell)
            nodes_per_edge
            for k=1:nodes_per_edge
                active_node((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=cell_j(i)*(nodes_per_edge-1)*nodes_per_row+cell_i(i)*(nodes_per_edge-1)+1+nodes_per_row*(k-1):cell_j(i)*(nodes_per_edge-1)*nodes_per_row+cell_i(i)*(nodes_per_edge-1)+nodes_per_edge+nodes_per_row*(k-1);
                local_node((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=nodes_per_edge*(k-1)+1:nodes_per_edge*k;
                active_cell_vector((i-1)*nodes_per_edge^2+nodes_per_edge*(k-1)+1:(i-1)*nodes_per_edge^2+k*nodes_per_edge)=active_cell(i);
            endfor
        endfor
        [active_node,local_node,active_cell_vector]=hf_remove_duplicates(active_node,local_node,active_cell_vector)
        for i=1:length(active_node)
            add=u(active_node(i))*hf_eval_poly_transformed(x,y,SF(local_node(i),:),mesh_size,vertex_matrix(cell_matrix(active_cell_vector(i),1),:));
            val+=add
        endfor
    endif


    
endfunction