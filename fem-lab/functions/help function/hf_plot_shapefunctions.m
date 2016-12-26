function hf_plot_shapefunctions()

  graphics_toolkit("gnuplot")
  addpath 'D:\Dokumente\GitHub\fem-lab\fem-lab\functions\shape function'
  addpath 'D:\Dokumente\GitHub\fem-lab\fem-lab\functions\mesh'
  
  for p=1:5
  
    polynomial_deg=p
    
    A = sf_generate(polynomial_deg);
    
    x=0:0.05:1;
    y=x;
    
    for k=1:rows(A);
    
      for i=1:length(x)
        for j=1:length(y)
          Z(i,j)=hf_eval_poly(x(i),y(j),A(k,:));
        endfor
      endfor
      
      h=figure(k)
      mesh(x,y,Z)
      axis ([0 1 0 1 0-(p-1)*0.25 1+(p-1)*0.25])
      xlabel ("x")
      ylabel ("y")
      zlabel ("p(x,y)")
      print(h,'-dpng',num2str(k),['D:\Dokumente\GitHub\fem-lab\fem-lab\functions\help function\shape functions of degree ' num2str(polynomial_deg) '\P' num2str(k)])
    
    endfor
  
  endfor
  
endfunction