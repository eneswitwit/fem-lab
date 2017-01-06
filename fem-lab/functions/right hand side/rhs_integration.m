function rhs = rhs_integration(mesh_size,polynomial_deg,cell,vertex,SF)

    cells_per_row=(1/mesh_size);
    number_of_nodes=((polynomial_deg*cells_per_row)+1)^2;
    
    for i=1:number_of_nodes
        RHS=0;
        [active_cell,active_sf]=node_node_to_cells(i,mesh_size,polynomial_deg);
        for k=1:rows(active_cell)
            active_cell_no=active_cell(k);
            active_vertex=vertex(active_cell_no,:);
            active_sf_no=active_sf(k);
            f = @(x,y) cos(mesh_size*x*pi+active_vertex(1))*cos(mesh_size*y*pi+active_vertex(2)).*hf_eval_poly(x,y,SF(active_sf_no,:));
            integral=dblquad(f,0,1,0,1,0.01);
            RHS=RHS+integral;
        endfor
        rhst(i)=RHS;
    endfor

    rhs=mesh_size^2*rhst';
    
endfunction