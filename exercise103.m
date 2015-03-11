%Exercise 10.3

%Values for Q
lambda_max = .1;
r_max = 2;
p_max = 3;
p_dot_max = 100000;

Q = [1/(lambda_max)^2 0 0 0;
     0 1/(r_max)^2 0 0;
     0 0 1/(p_max)^2 0;
     0 0 0 1/(p_dot_max)^2];
 
 Q = [10 0 0 0;
     0 .2 0 0;
     0 0 .1 0;
     0 0 0 .01];
 
 R = 1;
 
 %Finding K
 [K, S, E] = dlqr(A,B,Q,R);
 
 %u_k = u_k_ - K'*(x_k - x_k_);
 
