function m_plot_shape_functions()
    for p=1:5
        polynomial_deg=p;
        SF = sf_generate(polynomial_deg);
        x=0:0.05/p:1;
        y=x;
        [X,Y]=meshgrid(x,y);
        for k=1:rows(SF);
            Z=hf_eval_poly(X,Y,SF(k,:));
            h=figure(k);
            mesh(x,y,Z);
            axis ([0 1 0 1 0-(p-1)*0.25 1+(p-1)*0.25])
            xlabel ("x")
            ylabel ("y")
            zlabel ("p(x,y)")
            print(h,'-dpng',num2str(k),[pwd '\Shape Function Plots\Shape Functions of Degree ' num2str(polynomial_deg) '\P' num2str(k)])
         endfor
    endfor
endfunction