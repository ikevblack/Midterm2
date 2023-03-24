function [ x, niters ] = CG( A, b, x0 )                         %Algorithm from 8.3.5.2
r = b;                                                          %r^(0) := b
niters = 0;                                                     %k is niters, k := 0
xk = x0;                                                        %prepare first xk. algorithm says x^(0) = 0, but ignoring that since x0 is a parameter
while norm(r) > eps(1)*norm(b)                                  %r^(k) /= 0. Using stopping criteria from 8.3.6.1
    if niters > 1                                               %update x^(k)
        xk = x;
    end 
    if niters == 0                                              %if k == 0
        p = r;                                                  %p^(k) := r^(0)
    else
        gamma = -(transpose(p)*A*r)/(transpose(p)*A*p);         %gamma calculation
        p = r + gamma*p;                                        %p^(k) := r^(k) + gamma*p^(k-1)
    end
    alpha = (transpose(r)*r)/(transpose(p)*A*p);                %alpha calculatiom
    x = xk + alpha*p;                                           %x^(k+1) = x^(k) + alpha*p^(k)
    r = r - alpha*A*p;                                          %r^(k+1) = r^(k) - alpha*A*p^(k)
    niters = niters + 1;                                        %k := k + 1
end