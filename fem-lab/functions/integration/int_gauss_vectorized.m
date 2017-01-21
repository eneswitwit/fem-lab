function val = int_gauss_vectorized(x_1,w_1,x_2,w_2,f)         % Initialize matrices: one entry of X_1 and X_2 is exactly one possible combination of the quadrature points.     % computing f(X_1,X_2) will result in a matrix, with a value for each combination    % Same principle applies to gauss weights    [X_1,X_2]=meshgrid(x_1,x_2);    [W_1,W_2]=meshgrid(w_1,w_2);        % We initialize a matrix, with the values of f and the corresponding weights    value_matrix=f(X_1,X_2).*W_1.*W_1;    % The final value is obtained by adding all values of the matrix    val=sum(sum(value_matrix));       endfunction