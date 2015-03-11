%Exercise 2++

%Variables
T = 0.25;
q = 1;
N = 99;

%Continous model
A_c = [ 0   1   0           0;
        0   0   -K_2        0;
        0   0   0           1;
        0   0   -K_1*K_pp   -K_1*K_pd  ];
    
B_c = [ 0;
        0;
        0;
        K_1*K_pp];

%Discrete model    
A = A_c*T + eye(4);
B = B_c*T;

%Initial values
lambda_0 = pi;
r_0 = 0;
p_0 = 0;
p_dot_0 = 0;
x0 = [lambda_0, r_0, p_0, p_dot_0];

%QP matrices
G1 = [1, 0, 0, 0;
    0, 0, 0, 0;
    0, 0, q, 0;
    0, 0, 0, 0];
P1 = q;                                
G = 2*genq2(G1,P1,N,N,1);             
c = zeros(N*(size(A, 2) + size(B, 2)) ,1);                 

%Bounds
lowerBound = [
    repmat([-Inf; -Inf; -pi/6; -Inf], N, 1);
    repmat([-Inf], N, 1)];

upperBound = [
    repmat([Inf; Inf; pi/6; Inf], N, 1);
    repmat([Inf], N, 1)];

%Inequality constraints
A_ieq = zeros(N*(size(A, 2) + size(B, 2)), N*(size(A, 2) + size(B, 2)));
b_ieq = zeros(N*(size(A, 2) + size(B, 2)), 1);

%Equality constraints
A_eq = gena2(A, B, N, size(A, 2), size(B, 2));
b_eq = zeros(N*size(A, 2), 1);

%Adding final equaility constraint
A_eq = [A_eq; zeros(size(A, 1), (N-1)*(size(A, 2))) eye(size(A, 2)) zeros(size(B, 1), N*size(B, 2))];
b_eq = [b_eq; zeros(size(B, 1), size(B, 2))];

%Inital equality values
b_eq(1) = lambda_0;
b_eq(2) = r_0;
b_eq(3) = p_0;
b_eq(4) = p_dot_0;

%Solvin the QP-problem
z = quadprog(G, c, A_ieq, b_ieq, A_eq, b_eq, lowerBound, upperBound);

%free space reserved kjetil



%Plotting
delta_t = T;

mx = size(A,2);                        % Number of states (number of columns in A)
mu = size(B,2);                        % Number of inputs(number of columns in B)

M = N;

u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

length(u)

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

Antall = 5/delta_t;
Nuller = zeros(Antall,1);
Enere  = ones(Antall,1);

u   = [Nuller; u; Nuller];
x1  = [pi*Enere; x1; Nuller];
x2  = [Nuller; x2; Nuller];
x3  = [Nuller; x3; Nuller];
x4  = [Nuller; x4; Nuller];

t = 0:delta_t:delta_t*(length(u)-1);          

figure(2)
subplot(511)
stairs(t,u),grid
ylabel('u')
subplot(512)
plot(t,x1),grid
ylabel('lambda')
subplot(513)
plot(t,x2),grid
ylabel('r')
subplot(514)
plot(t,x3),grid
ylabel('p')
subplot(515)
plot(t,x4),grid
xlabel('tid (s)'),ylabel('pdot')
    
Nuller = zeros(Antall,4);
optimal_trajectory_data = [x1,x2,x3, x4];
optimal_trajectory.time = t;
optimal_trajectory.signals.values = optimal_trajectory_data;
optimal_trajectory.signals.dimension = 141;

p_c.time = t;
p_c.signals.values = u;
p_c.signals.dimension = 141;