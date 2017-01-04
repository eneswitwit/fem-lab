function test()

    graphics_toolkit("gnuplot")

    x=0:0.1:1
    y=0:0.1:1
    h=0.5
    
    for i = 1:length(x)
        for j=1:length(y)
            Z(i,j)=hf_eval_poly(x(i),y(j),[1 -1 -1 1]);
        endfor
    endfor
    
    mesh(x,y,Z)

end