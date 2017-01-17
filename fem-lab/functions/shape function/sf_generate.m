function SF = sf_generate(polynomial_deg)    % This script defines a matrix called SF saving in each row a shape function. The shape function is represented as coefficient vector.    % The shape functions are in lexicographical order. More over the shape functions are defined on the unit square.    % We want to use transformations to transform those shape functions to the shape functions of individual cells.        % ------ Mesh and cell properties -------    dimension_domain = 2;    k = 0;    nodes_counter =  (factorial(dimension_domain)/(factorial(k)*factorial(dimension_domain-k))) * polynomial_deg^k * (polynomial_deg+1)^(dimension_domain-k);    row_counter = polynomial_deg + 1 ;    nodes_per_row = nodes_counter / row_counter;    % Calculate mesh_size so we can generate a mesh. We can use our mesh_generate() function because we want our shape functions on the    % reference cell, which is a unit square. Our nodes should be equidistant. mesh_generate() perfectley does that.    mesh_size = 1/(nodes_per_row-1);    % Our nodes matrix is equal to the vertex matrix in this case    [node,cell] = mesh_generate(mesh_size);    % ----------------------------------------        % ---- Compute shape function matrix --------    for i=1:nodes_counter        % We need a matrix containing data for node i . I-th Shape function on node i equals to one, on other nodes 0.        data_set= zeros(length(node),3);        for j=1:length(node)            data_set(j,:) = [node(j,1),node(j,2), hf_kronecker_delta(j,i) ];        endfor                % Interpolate to get the i-th shape function.        SF(i,:) = hf_interpolation(data_set);    endfor    % -------------------------------------------    endfunction