function error_map()    % Add all subfolders to working directory    addpath(genpath([pwd '/functions']))        max_mesh_size_counter = 10;    max_pol_deg = 3;        error_map = zeros(max_mesh_size_counter, max_pol_deg);    for pol_deg = 1:max_pol_deg        for mesh_size_counter=1:max_mesh_size_counter            [error_L2,runtime] = main(1/(2*mesh_size_counter),pol_deg) ;            error_map(mesh_size_counter,pol_deg) = error_L2 ;         endfor    endfor    % Plot error map    for mesh_size_counter = 1:max_mesh_size_counter        x(mesh_size_counter) = 1/ (mesh_size_counter*2);    endfor    pol_deg = 1:max_pol_deg;    figure(4);    x    pol_deg    error_map    surf(x,pol_deg ,error_map');    set(gca,'zscale','log')    %axis([1 max_pol_deg 0 error_map(1,1) ]);    xlabel ("mesh size");    ylabel ("polynomial degree");    zlabel ("error");endfunction