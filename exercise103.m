%Exercise 10.3
Q = [10 0 0 0;
     0 .2 0 0;
     0 0 .1 0;
     0 0 0 .01];

 R = 1;
 
 %Finding K
 [K, S, E] = dlqr(A,B,Q,R);
 
 %u_k = u_k_ - K'*(x_k - x_k_);
 
