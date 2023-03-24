function [ x, niters ] = PCG( A, b, x0 )                                   %using algorithm from 8.3.6.1
L = ichol(sparse(A), struct('type','ict','droptol',1e-3,'michol','off'));  %using the following preconditioner
M = L*transpose(L);                                                        %establish M = L * L^T
xk = x0;                                                                   %prepare first xk
r = b;                                                                     %r^(0) := b
niters = 0;                                                                %niters is k. k := 0
while norm(r) > eps(1)*norm(b)                                             %r^(k) /= 0. Using stopping criteria from 8.3.6.1
    if niters > 1                                                          %update x^(k)
        xk = x;
    end 
    if niters == 0                                                         %if k == 0
        p = M\r;                                                           %p := (M^-1)*r^(0)
    else 
        gamma = (transpose(r)*(M\r))/(transpose(r_prev)*(M\r_prev));       %gamma calculation
        p = (M\r) + gamma*p;                                               %p^k := (M^-1) * r^(k) + gamma * p^(k-1)
    end
    q = A*p;                                                               %q := A * p^(k)
    alpha = (transpose(r)*(M\r))/(transpose(p)*q);                         %alpha calculation
    x = xk + alpha*p;                                                      %x^(k+1) := x^(k) + alpha*p^(k)
    r_prev = r;                                                            %preserve r^(k)
    r = r - alpha*q;                                                       %r^(k+1) := r^(k) - alpha*q^(k)
    niters = niters + 1;                                                   %k := k + 1
end