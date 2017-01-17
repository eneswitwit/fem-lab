function node = node_generate(polynomial_deg,mesh_size)    % This script computes a matrix saving all nodes in lexicographical order of the mesh globally.    % The first column represents the x value and the 2nd column the y value of a node.        % We start in the left corner of our unit square which is (0,0)    x=y=0;    node_nr=1;        % While we are not at the right corner of our unit square which is (1,1) we do the following:    while y <= 1 &&  x <= 1        % While we are not at the boundary in x direction we can walk to the right         while x<=1             node(node_nr,:)=[x,y];            x=x+mesh_size/polynomial_deg;            node_nr++;        endwhile        % We arrived at the boundary in x direction. Now we need to go one row up in y direction and go again to the first node from the left.        y=y+mesh_size/polynomial_deg;        x=0;            endwhile    endfunction