function rhs = rhs_integration(polynomial_deg,mesh_size,cell,vertix)

    cells_per_row=(1/mesh_size);
    number_of_nodes=((polynomial_deg*cells_per_row)+1)^2;
    
    for i=1:number_of_nodes
        RHS=0
        A=node_nodes_to_cells(i,mesh_size,polynomial_deg)
        active_cells=A(:,1)
        active_sf=A(:,2)    
        for k=1:rows(active_cells)
            active_cell_no=active_cells(k)
            active_sf_no=active_sf(k)
            f = @(x,y) cos(mesh_size*x*pi+vertix(active_cell,1)*cos(mesh_size*y*pi+vertix(active_cell,2).*hf_eval_poly(x,y,SF(active_sf_no))
            RHS=RHS+dblquad(f,0,1,0,1)
        endfor
        rhs(i)=RHS
    endfor

endfunction