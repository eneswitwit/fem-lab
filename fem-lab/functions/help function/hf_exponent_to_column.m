% Goal of this function is to put in the exponent of x and y and to get the index of the column in the coefficient vector belonging to these exponents.function index = hf_exponent_to_column(x_exponent,y_exponent)    max_combinations = factorial(x_exponent)*factorial(y_exponent);      data_set = zeros(max_combinations,3);    j=1;     % order of basis in the coefficient vector    level=0; % the level represents the structure of the ordering of the basisfunctions of the polynomial. It could be described as highest exponent allowed at a time.    while j<max_combinations        exponent_x=level;        exponent_y=0;        while(exponent_y < level)            data_set(j,:) = [exponent_x , exponent_y, j];            exponent_y++;            j++;        endwhile        exponent_x = 0;        exponent_y = level;        while(base_x < level)            data_set(j,:) = [exponent_x , exponent_y, j];            exponent_x++;            j++;        endwhile        data_set(j,:)=[level,level,j];        level++;        j++;    endwhile    % Check where x_exponent is in data_set    for i=1:length(data_set)        if(x_exponent == data_set(i,1))            % Check if y_exponent fits also.            if(y_exponent == data_set(i,2))                break;            endif        endif    endfor    index = data_set(i,3);endfunction