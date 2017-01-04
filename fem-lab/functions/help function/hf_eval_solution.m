function val=hf_eval_solution(x,y,u,mesh_size,polynomial_deg)

    cells_per_row=(1/mesh_size);
    nodes_per_row=(polynomial_deg*cells_per_row)+1;
    nodes_per_edge=polynomial_deg+1;
    
    if (floor(x*(nodes_per_row-1))==x*(nodes_per_row-1) && floor(y*(nodes_per_row-1))==y*(nodes_per_row-1))
        node_i=x*(nodes_per_row-1)
        node_j=y*(nodes_per_row-1)
        node_no=node_j*nodes_per_row+node_i+1
        %val=u(node_no)
    else
        if (floor(x*cells_per_row)==x*cells_per_row)        % x on edge
            cell_j=floor(y*cells_per_row)
            if(x==0)                                                    % x on lower boundary
                active_cell=(cell_j)*cells_per_row+1
            else                                                   % x not on lower boundary
                if(x==1)                                                    % x on upper boundary
                    active_cell=(cell_j+1)*cells_per_row
                else                                                   % x not on lower boundary -> on an interior edge
                    cell_i=x*cells_per_row-1
                    active_cell=[cell_j*cells_per_row+cell_i+1;cell_j*cells_per_row+cell_i+2]
                endif
            endif
        else                                                   % x not on edge
            if (floor(y*cells_per_row)==y*cells_per_row)            % y on edge
                cell_i=floor(x*cells_per_row)
                if (y==0)                               % y on lower boundary
                    active_cell=cell_i+1
                else                                                    % y not on lower boundary
                    if (y==1)                     % y on upper boundary
                        active_cell=cells_per_row*(cells_per_row-1)+cell_i+1
                    else                                                      % y not on lower boundary -> y on an interior edge
                        cell_j=y*cells_per_row-1
                        active_cell=[cell_j*cells_per_row+cell_i+1;(cell_j+1)*cells_per_row+cell_i+1]
                    endif
                endif
            else                                                    % y not on any edge -> (x,y) in the interior
                cell_i=floor(x*cells_per_row)
                cell_j=floor(y*cells_per_row)
                active_cell=cell_j*cells_per_row+cell_i+1
            endif
        endif
    endif
    
    
    
endfunction