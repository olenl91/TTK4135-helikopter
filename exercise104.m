clc

%exercise 10.4

%x = [lambda r p p_dot e e_dot]'
%u = [p_c e_c]'

% system parameters

q1 = 1;
q2 = q1;
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;
T = 0.25;
N = 40;
lambda_f = 0;

%initial values
lambda_0 = pi;
r_0 = 0;
p_0 = 0;
p_dot_0 = 0;
e_0 = 0;
e_dot_0 = 0;


A_c = [0 1 0 0 0 0;
        0 0 -K_2 0 0 0;
        0 0 0 1 0 0;
        0 0 -K_1*K_pp -K_1*K_pd 0 0;
        0 0 0 0 0 1;
        0 0 0 0 -K_3*K_ep -K_3*K_ed];

b_c = [0 0; 0 0; 0 0; K_1*K_pp 0; 0 0; 0 K_3*K_pp];

A = A_c*T + eye(6);
b = b_c*T;

% inequality constraints
%see mycon

% objective function
myf = @(lambda, p, e) ones(1, N)*(lambda.^2 + q1.*p.^2 + q2.*e.^2); % legg inn lambda_f ledd
f = @(x) myf(x(1: 6: 6*N), x(3 : 6: 6*N), x(5: 6: 6*N));

%Equality constraints
A_eq = gena2(A, b, N, size(A, 2), size(b, 2));
b_eq = zeros(N*size(A, 2), 1);

%Inital equality values
b_eq(1) = lambda_0;
b_eq(2) = r_0;
b_eq(3) = p_0;
b_eq(4) = p_dot_0;
b_eq(5) = e_0;
b_eq(6) = e_dot_0;

%x = [lambda r p p_dot e e_dot]'
%u = [p_c e_c]'

%Bounds
lowerBound = [
    repmat([-Inf; -Inf; -pi/6; -Inf; -Inf; -Inf], N, 1);
    repmat([-Inf; -Inf], N, 1)];

upperBound = [
    repmat([Inf; Inf; pi/6; Inf; Inf; Inf], N, 1);
    repmat([Inf; Inf], N, 1)];


z_0 = zeros((6+2)*40, 1);

OPTIONS = optimset('Algorithm','sqp');
z = FMINCON(f, z_0, [], [], A_eq, b_eq, lowerBound, upperBound);%, @mycon, OPTIONS);

