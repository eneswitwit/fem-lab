function correct = test_rhs_integration()    % Test whether rhs_integration delivers correct values by testing a special case.     % We will consider polynomial degree 1 and mesh size 0.5 and hardcode the solution for that case.    % correct will be a vector which will contain 0 and 1. If the entry is 0 at the i-th position, it means the vector delivered by     % rhs_integration is not correct in the i-th position, otherwise it is correct.        % Shape functions    % 1. Element    shape_function_1 = [1,-2,-2,4;0,2,0,-4;0,0,2,-4;0,0,0,4];    % 2. Element    shape_function_2 = [2,-2,-4,4;-1,2,2,-4;0,0,4,-4;0,0,-2,4];    % 3. Element    shape_function_3 = [2,-4,-2,4;0,4,0,-4;-1,2,2,-4;0,-2,0,4];    % 4. Element    shape_function_4 = [4,-4,-4,4;-2,4,2,-4;-2,2,4,-4;1,-2,-2,4];        % Initialize Gauss quadratures nodes and weights    [x_1,w_1] = int_gauss_weights(10,0,0.5);    [x_2,w_2] = int_gauss_weights(10,0.5,1);        % create right hand side    cell_matrix = mesh_cell(0.5,1);    rhs = zeros(9,1);         % -- First node --    shape_function_T1 = shape_function_1(1,:);    f_1 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T1);    %rhs(1) = dblquad(f_1,0,0.5,0,0.5);    rhs(1) = int_gauss_vectorized(x_1,w_1,x_1,w_1,f_1);        % -- Second node --    shape_function_T1 = shape_function_1(2,:);    shape_function_T2 = shape_function_2(1,:);    f_1 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T1);    f_2 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T2);    %rhs(2) = dblquad(f_1,0,0.5,0,0.5) + dblquad(f_2,0.5,1,0,0.5);    rhs(2) = int_gauss_vectorized(x_1,w_1,x_1,w_1,f_1) + int_gauss_vectorized(x_2,w_2,x_1,w_1,f_2);        % --Third node --    shape_function_T2 = shape_function_2(2,:);    f = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T2);    %rhs(3) = dblquad(f,0.5,1,0,0.5);    rhs(3) = int_gauss_vectorized(x_2,w_2,x_1,w_1,f);        % -- Fourth node --    shape_function_T1 = shape_function_1(3,:);    shape_function_T3 = shape_function_3(1,:);    f_1 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T1);    f_3 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T3);    %rhs(4) = dblquad(f_1,0,0.5,0,0.5) + dblquad(f_3,0,0.5,0.5,1);    rhs(4) = int_gauss_vectorized(x_1,w_1,x_1,w_1,f_1) + int_gauss_vectorized(x_1,w_1,x_2,w_2,f_3);        % -- Fifth node --    shape_function_T1 = shape_function_1(4,:);    shape_function_T2 = shape_function_2(3,:);    shape_function_T3 = shape_function_3(2,:);    shape_function_T4 = shape_function_4(1,:);    f_1 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T1);    f_2 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T2);    f_3 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T3);    f_4 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T4);    %rhs(5) = dblquad(f_1,0,0.5,0,0.5) + dblquad(f_3,0,0.5,0.5,1) + dblquad(f_2,0.5,1,0,0.5) + dblquad(f_4,0.5,1,0.5,1);    rhs(5) = int_gauss_vectorized(x_1,w_1,x_1,w_1,f_1) + int_gauss_vectorized(x_1,w_1,x_2,w_2,f_3) + int_gauss_vectorized(x_2,w_2,x_1,w_1,f_2) + int_gauss_vectorized(x_2,w_2,x_2,w_2,f_4);        % -- Sixth node --    shape_function_T2 = shape_function_2(4,:);    shape_function_T4 = shape_function_4(2,:);    f_2 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T2);    f_4 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T4);    %rhs(6) = dblquad(f_2,0.5,1,0,0.5) + dblquad(f_4,0.5,1,0.5,1);    rhs(6) = int_gauss_vectorized(x_2,w_2,x_1,w_1,f_2) + int_gauss_vectorized(x_2,w_2,x_2,w_2,f_4);    % -- Seventh node --    shape_function_T3 = shape_function_3(3,:);    f_3 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T3);    %rhs(7) = dblquad(f_3,0,0.5,0.5,1);    rhs(7) = int_gauss_vectorized(x_1,w_1,x_2,w_2,f_3);    % -- Eigth node --    shape_function_T3 = shape_function_3(4,:);    shape_function_T4 = shape_function_4(3,:);    f_3 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T3);    f_4 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T4);    %rhs(8) = dblquad(f_3,0,0.5,0.5,1) + dblquad(f_4,0.5,1,0.5,1);    rhs(8) = int_gauss_vectorized(x_1,w_1,x_2,w_2,f_3) + int_gauss_vectorized(x_2,w_2,x_2,w_2,f_4);    % -- Ninth node --    shape_function_T4 = shape_function_4(4,:);    f_4 = @(x,y) cos(pi.*x).*cos(pi.*y).*hf_eval_poly(x,y,shape_function_T4);    %rhs(9) = dblquad(f_4,0.5,1,0.5,1);    rhs(9) = int_gauss_vectorized(x_2,w_2,x_2,w_2,f_4);    % Generate values by rhs_integration, which we want to test by comparing it to the hardcoded solution    [Vertex,Cell] = mesh_generate(0.5);    SF = sf_generate(1);    rhs_generated_new = rhs_integration_new(Vertex,Cell,SF,@(x,y) cos(pi*x).*cos(pi*y))    rhs_generated = rhs_integration(Vertex,Cell,SF,@(x,y) cos(pi*x).*cos(pi*y))      % Test for correctness    correct = zeros(9,1);    for i=1:9        if( abs(rhs(i) - rhs_generated(i)) < 0.001 )            correct(i) = 1;        endif    endfor        % If everything is correct output just a 1 instead of the whole correct vector    test = 1;    for i=1:9        if( correct(i) == 0 )            test = 0;        endif    endfor    if( test == 1 )        correct = 1;    endif    endfunction