function val = error_L2(f,g)
    
    %this function will give us an estimation of the L2-error, by using Gauss-quadrature
    
    a=0;
    b=1;
    N=3;
    
    h=@(x,y) abs(f(x,y)-g(x,y)) 
    [x,w]=int_gauss_weights(N,a,b);
    val=int_gauss(x,w,x,w,h)
    
endfunction