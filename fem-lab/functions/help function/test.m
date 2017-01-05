function test()
    
    tic
    
    for k=1:5
    
    mesh_size=0.5;
    pol_deg=k;
    [vertex,cell]=mesh_generate(mesh_size);
    SF=sf_generate(pol_deg);
    cells_per_row=(1/mesh_size);
    nodes_per_row=(pol_deg*cells_per_row)+1;
       
    u=cos((1/nodes_per_row)*(1:nodes_per_row^2))  
    
    x=0:0.05:1;
    y=0:0.05:1;
    
    for i = 1:length(x)
        for j=1:length(y)
            Z(j,i)=hf_eval_solution(x(i),y(j),u,cell,vertex,pol_deg,SF);
        endfor
    endfor
    
    x
    y
    Z
    
    figure(k)
    mesh(x,y,Z)
    
    endfor
    
    toc
    
end