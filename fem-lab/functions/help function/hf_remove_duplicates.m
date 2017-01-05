function [active_node,local_node,active_cell_vector]=hf_remove_duplicates(active_node,local_node,active_cell_vector)

    for i=1:length(active_node)-1
        for j=i+1:length(active_node)
            if (active_node(j)==active_node(i))
                active_node(j)=0;
                local_node(j)=0;
                active_cell_vector(j)=0;
            endif
        end
    endfor
    
    Aux_matrix=unique([active_node',local_node',active_cell_vector'],"rows");
    if(Aux_matrix(1,1)==0)
        active_node=Aux_matrix(2:rows(Aux_matrix),1)';
        local_node=Aux_matrix(2:rows(Aux_matrix),2)';
        active_cell_vector=Aux_matrix(2:rows(Aux_matrix),3)';
    else
        active_node=Aux_matrix(1:rows(Aux_matrix),1)';
        local_node=Aux_matrix(1:rows(Aux_matrix),2)';
        active_cell_vector=Aux_matrix(1:rows(Aux_matrix),3)';
    endif
end