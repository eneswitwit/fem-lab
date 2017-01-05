function val = hf_eval_poly_transformed(x,y,coefficient_vector,scaling,displacement)
    j=1;        % order of basis in the coefficient vector
    level=0;    % the level represents the structure of the ordering of the basisfunctions of the polynomial. It could be described
    val=0;      % as the highest exponent allowed at a specific moment.
    while j < length(coefficient_vector)
        exponent_x=level;
        exponent_y=0;
        while(exponent_y < level)
            val= val + coefficient_vector(j)*(scaling^(-1)*(x-displacement(1))).^(exponent_x)*(scaling^(-1)*(y-displacement(2))).^(exponent_y);
            exponent_y++;
            j++;
        endwhile
        exponent_x = 0;
        exponent_y = level;
        while(exponent_x < level)
            val= val+ coefficient_vector(j)*(scaling^(-1)*(x-displacement(1))).^(exponent_x)*(scaling^(-1)*(y-displacement(2))).^(exponent_y);
            exponent_x++;
            j++;
        endwhile
        val= val+ coefficient_vector(j)*(scaling^(-1)*(x-displacement(1))).^(level)*(scaling^(-1)*(y-displacement(2))).^(level);
        j++;
        level++;
    endwhile
endfunction