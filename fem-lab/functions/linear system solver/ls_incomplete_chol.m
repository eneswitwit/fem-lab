function L=ls_incomplete_chol(A)

    n = length(A);
      for k=1:n
        A(k,k) = A(k,k)^(1/2);
        for i=(k+1):n
            if (A(i,k)!=0)
                A(i,k) = A(i,k)/A(k,k);            
            endif
        endfor
        for j=(k+1):n
            for i=j:n
                if (A(i,j)!=0)
                    A(i,j) = A(i,j)-A(i,k)*A(j,k);  
                endif
            endfor
        endfor
      endfor

        for i=1:n
            for j=i+1:n
                A(i,j) = 0;
            endfor
        endfor            

endfunction