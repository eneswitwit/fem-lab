function rhs = rhs_integration(mesh_size,polynomial_deg,cell,vertex,SF)
    % Initialize Gauss Quadratur
    degree_polynomial = sqrt(length(SF))-1;
    [sample_points,weights] = int_gauss_weights(polynomial_deg*10,0,1);
    
    
    cells_per_row=(1/mesh_size);
    number_of_nodes=((polynomial_deg*cells_per_row)+1)^2;
    
    for i=1:number_of_nodes
        RHS=0;
        [active_cell,active_sf]=node_node_to_cells(i,mesh_size,polynomial_deg);
        for k=1:rows(active_cell)
            active_cell_no=active_cell(k);
            active_vertex=vertex(cell(active_cell_no,1),:);
            active_sf_no=active_sf(k);
            f = @(x,y) cos(mesh_size*x*pi+active_vertex(1))*cos(mesh_size*y*pi+active_vertex(2)).*hf_eval_poly(x,y,SF(active_sf_no,:));
            %f = @(x,y) cos((1/mesh_size)*x*pi)*cos((1/mesh_size)*y*pi).*hf_eval_poly(x,y,SF(active_sf_no,:));
            % integral=dblquad(f,0,1,0,1);
            integral=int_gauss(sample_points,weights,f);
            RHS=RHS+integral;
        endfor
        rhst(i)=RHS;
    endfor

    rhs=mesh_size^2*rhst';
    
endfunction