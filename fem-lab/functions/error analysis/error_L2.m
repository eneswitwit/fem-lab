function val = error_L2(f,g,order)
    
    [x,w]=int_gauss_weights(order,0,1);
    h=@(x,y) (f(x,y)-g(x,y))^2;
    val=sqrt(int_gauss(x,w,x,w,h)) 
    
endfunction