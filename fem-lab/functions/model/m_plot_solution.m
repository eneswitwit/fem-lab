function m_plot_solution(solution_vector,polynomial_deg,mesh_size)    % We want to plot our approximated solution and the analytical solution.        % ---- Mesh and shape function properties -----    SF = sf_generate(polynomial_deg);    [vertex,cell] = mesh_generate(mesh_size);    % ---------------------------------------------        % Analytical solution    f = @(x,y) (1/(1+2*pi^2)).*cos(pi.*x).*cos(pi.*y);        % Right hand side of the initial problem    g = @(x,y) cos(pi.*x).*cos(pi.*y);        % Plotting properties    x=0:0.1:1;    y=x;        % Initialize matrices saving the value of the function for specific x and y values.    Z = zeros(length(x),length(y));    A = zeros(length(x),length(y));        % Calculate matrix saving the value of the approximated function for values x(i) and y(j)    for i=1:length(x)        for j=1:length(y)            Z(i,j)=hf_eval_solution(x(i),y(j),solution_vector,cell,vertex,polynomial_deg,SF);        endfor    endfor        % Calculate matrices saving the value of the for analytical solution and right hand side of initial problem for values x(i) and y(j)    for i=1:length(x)        for j=1:length(y)            % Right hand side of initial problem            A(j,i)= g(x(i),y(j));            % Analytical solution            B(j,i)=f(x(i),y(j));        endfor    endfor        % Plot analytical solution    figure(1);    mesh(x,y,B);    xlabel ("x");    ylabel ("y");    zlabel ("u(x,y)");        % Plot right hand side of initial problem    figure(2);    mesh(x,y,A);    xlabel ("x");    ylabel ("y");    zlabel ("f(x,y)");        % Plot approximated solution by FEM    figure(3);    mesh(x,y,Z);    xlabel ("x");    ylabel ("y");    zlabel ("u_h(x,y)");    endfunction