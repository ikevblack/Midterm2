function [ x, niters ] = Method_of_Steepest_Descent( A, b, x0 ) %referencing Figure 8.2.4.1
  niters = 0;                                                   %niters is k, k := 0
  r = b - A*x0;                                                 %r0 := b - A*x0 
  xk = x0;                                                      %prepare first xk
  while norm(r) > eps(1)*(norm(b))                              %r^k /= 0. Using stopping criteria from 8.3.6.1
    if niters > 1                                               %update x^k
        xk = x;
    end 
    p = r;                                                      %p^k := r^k
    q = A * p;                                                  %q^k := A*p^k
    alpha = (transpose(p) * r)/(transpose(p) * q);              %alpha value
    x = xk + alpha*p;                                           %x^(k+1) := x^k + alpha*p^k
    r = r - alpha*q;                                            %r^(k+1) := r^k - alpha*q^k
    niters = niters + 1;                                        %k := k + 1
  end