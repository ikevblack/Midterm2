function [ x, niters ] = Method_of_Steepest_Descent_ichol( A, b, x0 )      %referencing right side equation of 8.2.5.1
L = ichol(sparse(A), struct('type','ict','droptol',6e-4,'michol','off'));  %using the following preconditioner
M = L*transpose(L);                                                        %establish M = LL^T
r = b - A*x0;                                                              %r := b - A*x0
xk = x0;                                                                   %prepare first xk
niters = 0;                                                                %niters is k, k := 0
while norm(r) > eps(1)*(norm(b))                                           %r^k /= 0. Using stopping criteria from 8.3.6.1
    if niters > 1                                                          %update x^k
        xk = x;
    end 
    p = M\r;                                                               %p^k := (M^-1)*r^k
    q = A*p;                                                               %q^k := A*p^k
    alpha = (transpose(p)*r)/(transpose(p)*q);                             %alpha calculation
    x = xk + alpha*p;                                                      %x^(k+1) := x^k + alpha*p^k
    r = r - alpha*q;                                                       %r^(K+1) := r^k - alpha*q^k
    niters = niters + 1;                                                   %k := k + 1
end